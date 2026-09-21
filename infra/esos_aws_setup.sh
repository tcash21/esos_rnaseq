#!/usr/bin/env bash
# =============================================================================
# ESOS genomics — AWS setup (run from your Mac, with AWS CLI v2 configured)
#   Creates: private encrypted S3 bucket, IAM role (S3 + SSM), no-inbound
#   security group, and one EC2 analysis instance with a 1 TB data volume.
#   Also adds an idle alarm that auto-STOPS the instance after ~1 hr of <3% CPU.
#
# Prereqs on the Mac:
#   brew install awscli && aws configure
#   brew install --cask session-manager-plugin     # to connect without SSH
#
# Put esos_bootstrap.sh in the same folder, then:  bash esos_aws_setup.sh
# =============================================================================
set -euo pipefail

REGION=us-east-1
# Reuse an existing bucket on re-runs:  BUCKET=esos-sema4-xxxxxxxx bash esos_aws_setup.sh
BUCKET="${BUCKET:-$(aws s3api list-buckets --query "Buckets[?starts_with(Name,'esos-sema4-')].Name | [0]" --output text 2>/dev/null)}"
[[ -z "$BUCKET" || "$BUCKET" == "None" ]] && BUCKET="esos-sema4-$(openssl rand -hex 4)"   # globally unique
INSTANCE_TYPE=m7i.4xlarge                    # 16 vCPU / 64 GB RAM
DATA_GB=1000                                 # /data volume (gp3)
ROLE=esos-ec2-role
PROFILE=esos-ec2-profile
TAG=esos-analysis

echo ">> Bucket: $BUCKET"

# ---------- 1. Private, encrypted S3 bucket ----------------------------------
if aws s3api head-bucket --bucket "$BUCKET" 2>/dev/null; then
  echo ">> Bucket exists, reusing"
else
  aws s3api create-bucket --bucket "$BUCKET" --region "$REGION"
fi
aws s3api put-public-access-block --bucket "$BUCKET" \
  --public-access-block-configuration \
  BlockPublicAcls=true,IgnorePublicAcls=true,BlockPublicPolicy=true,RestrictPublicBuckets=true
aws s3api put-bucket-encryption --bucket "$BUCKET" \
  --server-side-encryption-configuration \
  '{"Rules":[{"ApplyServerSideEncryptionByDefault":{"SSEAlgorithm":"AES256"}}]}'
# Everything moves to Intelligent-Tiering (cheaper when BAMs sit untouched)
aws s3api put-bucket-lifecycle-configuration --bucket "$BUCKET" \
  --lifecycle-configuration '{"Rules":[{"ID":"it","Status":"Enabled","Filter":{},
  "Transitions":[{"Days":0,"StorageClass":"INTELLIGENT_TIERING"}]}]}'

# ---------- 2. IAM role: this bucket only + SSM Session Manager --------------
cat > /tmp/esos-trust.json <<'JSON'
{"Version":"2012-10-17","Statement":[{"Effect":"Allow",
 "Principal":{"Service":"ec2.amazonaws.com"},"Action":"sts:AssumeRole"}]}
JSON
cat > /tmp/esos-s3.json <<JSON
{"Version":"2012-10-17","Statement":[
 {"Effect":"Allow","Action":["s3:ListBucket"],"Resource":"arn:aws:s3:::$BUCKET"},
 {"Effect":"Allow","Action":["s3:GetObject","s3:PutObject","s3:DeleteObject"],
  "Resource":"arn:aws:s3:::$BUCKET/*"}]}
JSON
aws iam get-role --role-name "$ROLE" >/dev/null 2>&1 || \
  aws iam create-role --role-name "$ROLE" \
    --assume-role-policy-document file:///tmp/esos-trust.json >/dev/null
aws iam put-role-policy --role-name "$ROLE" --policy-name esos-s3 \
  --policy-document file:///tmp/esos-s3.json          # overwrites -> points at $BUCKET
aws iam attach-role-policy --role-name "$ROLE" \
  --policy-arn arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore
aws iam get-instance-profile --instance-profile-name "$PROFILE" >/dev/null 2>&1 || \
  aws iam create-instance-profile --instance-profile-name "$PROFILE" >/dev/null
IN_PROFILE=$(aws iam get-instance-profile --instance-profile-name "$PROFILE" \
  --query "InstanceProfile.Roles[?RoleName=='$ROLE'] | length(@)" --output text)
[[ "$IN_PROFILE" == "0" ]] && \
  aws iam add-role-to-instance-profile --instance-profile-name "$PROFILE" --role-name "$ROLE"
echo ">> waiting for IAM to propagate..."; sleep 15

# ---------- 3. Network: default VPC + no-inbound security group -------------
VPC_ID=$(aws ec2 describe-vpcs --region "$REGION" \
  --filters Name=isDefault,Values=true --query 'Vpcs[0].VpcId' --output text)
if [[ -z "$VPC_ID" || "$VPC_ID" == "None" ]]; then
  echo ">> No default VPC in $REGION — creating one"
  VPC_ID=$(aws ec2 create-default-vpc --region "$REGION" --query 'Vpc.VpcId' --output text)
  sleep 10
fi
echo ">> VPC: $VPC_ID"
SG_ID=$(aws ec2 describe-security-groups --region "$REGION" \
  --filters Name=group-name,Values=esos-no-inbound Name=vpc-id,Values="$VPC_ID" \
  --query 'SecurityGroups[0].GroupId' --output text)
if [[ -z "$SG_ID" || "$SG_ID" == "None" ]]; then
  SG_ID=$(aws ec2 create-security-group --region "$REGION" \
    --group-name esos-no-inbound --description "No inbound; SSM only" \
    --vpc-id "$VPC_ID" --query GroupId --output text)
fi

# ---------- 4. Launch the instance -------------------------------------------
IID=$(aws ec2 run-instances --region "$REGION" \
  --image-id resolve:ssm:/aws/service/canonical/ubuntu/server/24.04/stable/current/amd64/hvm/ebs-gp3/ami-id \
  --instance-type "$INSTANCE_TYPE" \
  --iam-instance-profile Name="$PROFILE" \
  --security-group-ids "$SG_ID" \
  --metadata-options HttpTokens=required \
  --block-device-mappings "[
    {\"DeviceName\":\"/dev/sda1\",\"Ebs\":{\"VolumeSize\":100,\"VolumeType\":\"gp3\",\"Encrypted\":true,\"DeleteOnTermination\":true}},
    {\"DeviceName\":\"/dev/sdf\",\"Ebs\":{\"VolumeSize\":$DATA_GB,\"VolumeType\":\"gp3\",\"Iops\":3000,\"Throughput\":250,\"Encrypted\":true,\"DeleteOnTermination\":true}}]" \
  --user-data file://esos_bootstrap.sh \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$TAG}]" \
  --query 'Instances[0].InstanceId' --output text)
echo ">> Instance: $IID"

# ---------- 5. Cost guard: auto-stop after ~60 min idle ----------------------
aws cloudwatch put-metric-alarm --region "$REGION" \
  --alarm-name "esos-idle-stop-$IID" \
  --namespace AWS/EC2 --metric-name CPUUtilization \
  --dimensions Name=InstanceId,Value="$IID" \
  --statistic Average --period 300 --evaluation-periods 12 \
  --threshold 3 --comparison-operator LessThanThreshold \
  --alarm-actions "arn:aws:automate:${REGION}:ec2:stop"

cat <<EOF

Done.
  Bucket:   s3://$BUCKET
  Instance: $IID  ($INSTANCE_TYPE)

Setup takes ~10 min after boot (check: sudo tail -f /var/log/esos-bootstrap.log)

Connect:
  aws ssm start-session --region $REGION --target $IID
  sudo su - ubuntu

Stop when not using (EBS still billed):
  aws ec2 stop-instances --region $REGION --instance-ids $IID
Terminate when done (after results are in S3):
  aws ec2 terminate-instances --region $REGION --instance-ids $IID
EOF
