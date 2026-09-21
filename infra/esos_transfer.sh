#!/usr/bin/env bash
# =============================================================================
# esos_transfer.sh — Dropbox -> EC2 transfer of the Sema4 delivery
#
# BEFORE STARTING (from your Mac): pause the idle auto-stop alarm — a download
# uses almost no CPU, so the alarm would otherwise stop the instance mid-transfer:
#   aws cloudwatch disable-alarm-actions --region us-east-1 \
#     --alarm-names $(aws cloudwatch describe-alarms --region us-east-1 \
#       --alarm-name-prefix esos-idle-stop --query 'MetricAlarms[].AlarmName' --output text)
# (re-enable afterwards with enable-alarm-actions, same arguments)
#
# Run ON THE EC2 INSTANCE, inside tmux so it survives disconnects:
#   tmux new -s dl
#   bash esos_transfer.sh "https://www.dropbox.com/scl/fi/....?rlkey=...&dl=0" [--extract] [--s3 s3://BUCKET/sema4/]
#
#   --extract     also unpack into /data/sema4 after verifying
#   --s3 URI      also sync the extracted files to S3 (implies --extract)
#
# What it does:
#   1. forces the Dropbox link to direct-download (dl=1)
#   2. checks free disk space
#   3. downloads with automatic resume/retry (safe to re-run if interrupted)
#   4. verifies: size matches Dropbox, file isn't an HTML error page,
#      and the full archive reads cleanly (gzip CRC + tar structure)
#   5. saves a listing + a summary of what's inside, by file type
#   6. optionally extracts and syncs to S3
#
# Output (in /data): sema4.tar, sema4_listing.txt, sema4_summary.txt, transfer.log
# =============================================================================
set -euo pipefail

URL="${1:-}"; shift || true
EXTRACT=0; S3_URI=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --extract) EXTRACT=1 ;;
    --s3) S3_URI="$2"; EXTRACT=1; shift ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac; shift
done
[[ -z "$URL" ]] && { echo "Usage: $0 <dropbox-share-link> [--extract] [--s3 s3://bucket/prefix/]"; exit 1; }

DATA=/data
OUT="$DATA/sema4.tar"
EXTRACT_DIR="$DATA/sema4"
LOG="$DATA/transfer.log"
MIN_FREE_GB=250          # ~115 GB archive + ~115 GB extracted + headroom
MAX_ATTEMPTS=20

exec > >(tee -a "$LOG") 2>&1
log() { echo "[$(date '+%F %T')] $*"; }

# ---------- 1. Force direct download ----------
URL="${URL//dl=0/dl=1}"
[[ "$URL" == *"dl=1"* ]] || URL="${URL}$([[ "$URL" == *\?* ]] && echo '&' || echo '?')dl=1"
log "Starting transfer"

# ---------- 2. Disk space ----------
FREE_GB=$(df -BG --output=avail "$DATA" | tail -1 | tr -dc '0-9')
log "Free space on $DATA: ${FREE_GB} GB"
if (( FREE_GB < MIN_FREE_GB )) && [[ ! -f "$OUT" ]]; then
  log "ERROR: need at least ${MIN_FREE_GB} GB free"; exit 1
fi

# ---------- 3. Expected size from Dropbox ----------
REMOTE_BYTES=$(curl -sIL "$URL" | tr -d '\r' | awk 'tolower($1)=="content-length:"{v=$2} END{print v+0}')
if (( REMOTE_BYTES > 1000000 )); then
  log "Remote size: $(numfmt --to=iec "$REMOTE_BYTES") ($REMOTE_BYTES bytes)"
else
  log "WARN: Dropbox didn't report a size; will verify by reading the archive instead"
  REMOTE_BYTES=0
fi

local_bytes() { [[ -f "$OUT" ]] && stat -c %s "$OUT" || echo 0; }

# ---------- 4. Download with resume ----------
attempt=0
while :; do
  have=$(local_bytes)
  if (( REMOTE_BYTES > 0 && have == REMOTE_BYTES )); then
    log "Download complete ($(numfmt --to=iec "$have"))"; break
  fi
  (( attempt++ )) || true
  (( attempt > MAX_ATTEMPTS )) && { log "ERROR: gave up after $MAX_ATTEMPTS attempts"; exit 1; }
  log "Attempt $attempt — resuming from $(numfmt --to=iec "$have")"
  if curl -L --fail -C - --retry 5 --retry-delay 15 --connect-timeout 30 \
          --speed-limit 1048576 --speed-time 120 \
          -o "$OUT" "$URL"; then
    # If Dropbox gave no size, a clean curl exit is our completion signal
    (( REMOTE_BYTES == 0 )) && { log "Download finished"; break; }
  else
    log "curl exited with $? — retrying in 30s"; sleep 30
  fi
done

# ---------- 5. Verify ----------
FTYPE=$(file -b "$OUT")
log "File type: $FTYPE"
if [[ "$FTYPE" == *HTML* || $(local_bytes) -lt 1000000 ]]; then
  log "ERROR: got an HTML page, not the data. The link is probably not public,"
  log "       expired, or hit Dropbox's download limit. Check the share settings."
  exit 1
fi

TAR_Z=""
[[ "$FTYPE" == *gzip* ]] && TAR_Z="--use-compress-program=pigz"

log "Reading the full archive to verify integrity + build listing (~20-40 min)..."
if ! tar $TAR_Z -tvf "$OUT" > "$DATA/sema4_listing.txt"; then
  log "ERROR: archive is corrupt or truncated. Delete $OUT and re-run."; exit 1
fi
log "Archive OK — $(wc -l < "$DATA/sema4_listing.txt") entries"

# Summary by file type (handles .fastq.gz, .vcf.gz, .bam, etc.)
{
  echo "=== Contents by file type ==="
  awk '$1 !~ /^d/ {
         name=$6; for(i=7;i<=NF;i++) name=name" "$i
         n=split(name, p, "/"); f=p[n]
         ext=f; sub(/^[^.]*/, "", ext); if (ext=="") ext="(none)"
         cnt[ext]++; sz[ext]+=$3 }
       END { for (e in cnt) printf "%-28s %6d files  %10.2f GB\n", e, cnt[e], sz[e]/1e9 }' \
    "$DATA/sema4_listing.txt" | sort -k4 -nr
  echo
  echo "=== Largest 25 files ==="
  sort -k3 -nr "$DATA/sema4_listing.txt" | head -25 | awk '{printf "%10.2f GB  %s\n", $3/1e9, $NF}'
  echo
  echo "=== Likely RNA files ==="
  grep -Ei 'rna|wts|transcript' "$DATA/sema4_listing.txt" | awk '{print $NF}' || echo "(none matched)"
} | tee "$DATA/sema4_summary.txt"

# ---------- 6. Optional extract + S3 ----------
if (( EXTRACT )); then
  mkdir -p "$EXTRACT_DIR"
  log "Extracting to $EXTRACT_DIR ..."
  tar $TAR_Z -xf "$OUT" -C "$EXTRACT_DIR"
  log "Extracted: $(du -sh "$EXTRACT_DIR" | cut -f1)"
fi

if [[ -n "$S3_URI" ]]; then
  log "Syncing to $S3_URI ..."
  aws s3 sync "$EXTRACT_DIR" "$S3_URI" --only-show-errors
  aws s3 cp "$DATA/sema4_listing.txt" "${S3_URI%/}/_sema4_listing.txt" --only-show-errors
  log "S3 sync complete"
fi

log "DONE. Summary: $DATA/sema4_summary.txt"
log "Reminder: turn off the Dropbox share link now that the transfer is done."
