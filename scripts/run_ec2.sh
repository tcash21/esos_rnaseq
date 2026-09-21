#!/usr/bin/env bash
# =============================================================================
# run_ec2.sh — run one numbered analysis step on the EC2 instance, from the Mac.
#
#   bash scripts/run_ec2.sh 01          # runs scripts/01_*.sh then scripts/01_*.R on EC2
#   bash scripts/run_ec2.sh 01 --dry    # show what would run
#
# Why EC2: the Sema4 data already lives there (/data/sema4), bcftools is
# installed, and it has unrestricted network to ClinVar / gnomAD / UCSC.
# Flow: start instance if stopped -> push scripts/ to S3 -> remote pulls them,
# runs the step with PROJ=/data/proj (data/sema4 symlinked) -> results synced
# to S3 -> pulled back into ./results/ here. Non-interactive (SSM RunCommand),
# so it works even if your terminal disconnects.
# =============================================================================
set -euo pipefail

STEP="${1:-}"; [[ -n "$STEP" ]] || { echo "usage: $0 <step-number, e.g. 01> [--dry]"; exit 1; }
DRY="${2:-}"
REGION="${ESOS_REGION:-us-east-1}"
IID="${ESOS_ID:-i-09847e2ec5b85a3bc}"
BUCKET="${ESOS_BUCKET:-s3://esos-sema4-576d839b}"
PROJ_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJ_DIR"

SH=$(ls scripts/${STEP}_*.sh 2>/dev/null | head -1 || true)
RS=$(ls scripts/${STEP}_*.R  2>/dev/null | head -1 || true)
[[ -n "$SH$RS" ]] || { echo "no scripts/${STEP}_* found"; exit 1; }
echo ">> step $STEP: ${SH:-} ${RS:-}"

# Remote script (runs as root via SSM; conda envs live under /home/ubuntu)
read -r -d '' REMOTE <<EOF || true
set -euo pipefail
export PATH=/home/ubuntu/miniforge3/envs/bio/bin:/home/ubuntu/miniforge3/bin:/usr/local/bin:\$PATH
export HOME=/root TMPDIR=/data/tmp
mkdir -p /data/proj/data /data/proj/results/logs /data/tmp && cd /data/proj
[[ -e data/sema4 ]] || ln -s /data/sema4 data/sema4
aws s3 sync ${BUCKET}/proj/scripts/ scripts/ --only-show-errors
command -v Rscript >/dev/null || { apt-get update -qq && apt-get install -y -qq r-base-core r-cran-data.table >/dev/null; }
export PROJ=/data/proj
LOG=results/logs/step_${STEP}_\$(date +%Y%m%d_%H%M%S).log
{
  $( [[ -n "$SH" ]] && echo "echo '### bash $SH'; bash $SH" )
  $( [[ -n "$RS" ]] && echo "echo '### Rscript $RS'; Rscript $RS" )
} 2>&1 | tee "\$LOG" || RC=\$?
aws s3 sync results/ ${BUCKET}/proj/results/ --only-show-errors
echo "REMOTE DONE rc=\${RC:-0}"; exit \${RC:-0}
EOF

if [[ "$DRY" == "--dry" ]]; then echo "$REMOTE"; exit 0; fi

# ---- start instance if needed ----
state=$(aws ec2 describe-instances --region "$REGION" --instance-ids "$IID" \
  --query 'Reservations[0].Instances[0].State.Name' --output text)
if [[ "$state" != "running" ]]; then
  echo ">> instance is $state — starting"
  aws ec2 start-instances --region "$REGION" --instance-ids "$IID" >/dev/null
  aws ec2 wait instance-running --region "$REGION" --instance-ids "$IID"
  echo ">> waiting for SSM agent"; sleep 60
fi

# ---- push scripts, launch ----
aws s3 sync scripts/ "$BUCKET/proj/scripts/" --only-show-errors --exclude "*.log"
CMD_ID=$(aws ssm send-command --region "$REGION" --instance-ids "$IID" \
  --document-name AWS-RunShellScript --timeout-seconds 7200 \
  --parameters "$(python3 -c 'import json,sys; print(json.dumps({"commands":[sys.stdin.read()],"executionTimeout":["7200"]}))' <<< "$REMOTE")" \
  --query 'Command.CommandId' --output text)
echo ">> running on EC2 (command $CMD_ID) — this can take 10–30 min; safe to Ctrl-C, results still sync to S3"

# ---- poll ----
while :; do
  sleep 20
  st=$(aws ssm get-command-invocation --region "$REGION" --command-id "$CMD_ID" --instance-id "$IID" \
        --query 'Status' --output text 2>/dev/null || echo Pending)
  printf '   %s  %s\n' "$(date +%H:%M:%S)" "$st"
  case "$st" in Success|Failed|Cancelled|TimedOut) break ;; esac
done
aws ssm get-command-invocation --region "$REGION" --command-id "$CMD_ID" --instance-id "$IID" \
  --query '[StandardOutputContent,StandardErrorContent]' --output text | tail -60

# ---- pull results back ----
mkdir -p results
aws s3 sync "$BUCKET/proj/results/" results/ --only-show-errors
echo ">> results synced to $PROJ_DIR/results/   (status: $st)"
[[ "$st" == "Success" ]]
