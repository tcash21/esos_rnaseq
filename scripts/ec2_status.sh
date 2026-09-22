#!/usr/bin/env bash
# ec2_status.sh — one-shot progress check of whatever step is running on the EC2 instance, from the Mac.
#   bash scripts/ec2_status.sh          # instance state, latest step log tail, STAR progress, load
set -euo pipefail
REGION="${ESOS_REGION:-us-east-1}"; IID="${ESOS_ID:-i-09847e2ec5b85a3bc}"
state=$(aws ec2 describe-instances --region "$REGION" --instance-ids "$IID" --query 'Reservations[0].Instances[0].State.Name' --output text)
echo "instance: $state"
[[ "$state" == "running" ]] || exit 0
alarm=$(aws cloudwatch describe-alarms --region "$REGION" --alarm-names "esos-idle-stop-$IID" --query 'MetricAlarms[0].[ActionsEnabled,StateValue]' --output text)
echo "idle-stop alarm (actions enabled, state): $alarm"
CID=$(aws ssm send-command --region "$REGION" --instance-ids "$IID" --document-name AWS-RunShellScript \
  --parameters '{"commands":["#!/bin/bash","uptime | sed \"s/.*up/up/\"","L=$(ls -t /data/proj/results/logs/step_*.log 2>/dev/null | head -1); echo \"== $L\"; tail -n 6 \"$L\" | cut -c1-160","P=$(ls -t /data/proj/results/*/requant/*.Log.progress.out /data/proj/results/*/*.Log.progress.out 2>/dev/null | head -1); [[ -n \"$P\" ]] && { echo \"== STAR progress\"; tail -n 3 \"$P\"; }","pgrep -a STAR | head -1 | cut -c1-80; df -h /data | tail -1"]}' \
  --query 'Command.CommandId' --output text)
for i in $(seq 1 15); do
  sleep 2
  st=$(aws ssm get-command-invocation --region "$REGION" --command-id "$CID" --instance-id "$IID" --query 'Status' --output text 2>/dev/null || echo Pending)
  [[ "$st" == "Success" || "$st" == "Failed" ]] && break
done
aws ssm get-command-invocation --region "$REGION" --command-id "$CID" --instance-id "$IID" --query 'StandardOutputContent' --output text
