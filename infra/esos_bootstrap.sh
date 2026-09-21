#!/bin/bash
# =============================================================================
# EC2 user-data — runs once as root on first boot (Ubuntu 24.04).
# Log: /var/log/esos-bootstrap.log     Marker when finished: /data/BOOTSTRAP_DONE
# =============================================================================
set -euxo pipefail
exec > >(tee -a /var/log/esos-bootstrap.log) 2>&1
export DEBIAN_FRONTEND=noninteractive

# ---------- System packages ----------
apt-get update -y
apt-get install -y build-essential curl wget unzip pigz pv htop tmux git jq \
  openjdk-17-jre-headless docker.io
systemctl enable --now docker
usermod -aG docker ubuntu
snap install aws-cli --classic

# ---------- Format + mount the 1 TB data volume at /data ----------
ROOT_DISK=$(lsblk -no PKNAME "$(findmnt -no SOURCE /)")
DATA_DEV=$(lsblk -dpno NAME,TYPE | awk '$2=="disk"{print $1}' | grep -v "/dev/${ROOT_DISK}$" | head -1)
blkid "$DATA_DEV" || mkfs.ext4 -m 0 "$DATA_DEV"
mkdir -p /data
UUID=$(blkid -s UUID -o value "$DATA_DEV")
grep -q "$UUID" /etc/fstab || echo "UUID=$UUID /data ext4 defaults,nofail 0 2" >> /etc/fstab
mount -a
mkdir -p /data/{sema4,refs,work,results,tmp}
chown -R ubuntu:ubuntu /data

# ---------- Nextflow (for nf-core pipelines, runs with -profile docker) ----------
cd /tmp && curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/ && chmod 755 /usr/local/bin/nextflow

# ---------- Miniforge + bioinformatics envs (as ubuntu) ----------
UH=/home/ubuntu
sudo -u ubuntu bash -c "
  set -e
  cd $UH
  curl -sL -o mf.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash mf.sh -b -p $UH/miniforge3 && rm mf.sh
  $UH/miniforge3/bin/conda init bash
  $UH/miniforge3/bin/conda config --add channels bioconda
  $UH/miniforge3/bin/conda config --add channels conda-forge
  $UH/miniforge3/bin/conda config --set channel_priority strict
"
# Core tools: BAM/VCF handling, RNA quantification, alignment, fusions, coverage
sudo -u ubuntu $UH/miniforge3/bin/mamba create -y -n bio \
  samtools bcftools htslib salmon star arriba mosdepth \
  || echo "WARN: bio env failed — rerun manually"
# Allele-specific copy number (FACETS) — separate env to keep solves simple
sudo -u ubuntu $UH/miniforge3/bin/mamba create -y -n facets \
  snp-pileup r-facets \
  || echo "WARN: facets env failed — rerun manually"

# ---------- Shell defaults ----------
cat >> $UH/.bashrc <<'RC'
export TMPDIR=/data/tmp
export NXF_HOME=/data/work/.nextflow
cd /data
RC
chown ubuntu:ubuntu $UH/.bashrc

touch /data/BOOTSTRAP_DONE
echo "BOOTSTRAP COMPLETE"

# =============================================================================
# RUNBOOK (run manually after connecting; not executed here)
# -----------------------------------------------------------------------------
# 1) Pull the tar from Dropbox (share link, change dl=0 -> dl=1). Run in tmux:
#      tmux new -s dl
#      cd /data && curl -L -o sema4.tar "https://www.dropbox.com/...&dl=1"
#
# 2) See what's inside, then extract (tar auto-detects gzip):
#      tar -tvf sema4.tar | tee sema4_listing.txt | less
#      tar -xvf sema4.tar -C /data/sema4
#
# 3) Check the reference build / contig naming (report says GRCh37/hg19):
#      samtools view -H /data/sema4/<tumor>.bam | grep '^@SQ' | head
#
# 4) Back up extracted files to S3, then delete the tar locally:
#      aws s3 sync /data/sema4 s3://<BUCKET>/sema4/
#      rm /data/sema4.tar
#
# 5) Push results as you go:
#      aws s3 sync /data/results s3://<BUCKET>/results/
# =============================================================================
