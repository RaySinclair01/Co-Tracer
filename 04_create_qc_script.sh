#!/bin/bash
# 文件名: create_qc_script.sh
# 用途: 创建质控脚本并执行

# 确保scripts目录存在
mkdir -p ./scripts

cat > ./scripts/run_qc.sh << 'EOF'
#!/bin/bash

# 定义变量
R1="./raw_data/TZ_16-HDS1-40.R1.raw.fastq.gz"
R2="./raw_data/TZ_16-HDS1-40.R2.raw.fastq.gz"
QC_R1="./qc_data/L1HII0100002-TS_T1.R1.clean.fastq.gz"
QC_R2="./qc_data/L1HII0100002-TS_T1.R2.clean.fastq.gz"

# 创建必要的目录
mkdir -p ./qc_data
mkdir -p ./logs

# 检查原始数据是否存在
if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "错误：原始数据文件不存在，请检查 raw_data 目录"
    exit 1
fi

# 双端质控 - 使用fastp
echo "开始双端数据质控..."
fastp -i ${R1} -I ${R2} \
      -o ${QC_R1} -O ${QC_R2} \
      -h ./qc_data/fastp_report.html -j ./qc_data/fastp.json \
      --detect_adapter_for_pe \
      --cut_front --cut_tail \
      --cut_mean_quality 20 \
      --qualified_quality_phred 20 \
      --length_required 50 \
      --thread 38 \
      --compression 4 \
      --report_title "L1HII0100002-TS_T1 QC Report" \
      > ./logs/fastp.log 2>&1

# 检查质控结果
if [ -f "$QC_R1" ] && [ -f "$QC_R2" ]; then
    echo "质控成功完成。检查结果文件:"
    ls -lh ./qc_data/
else
    echo "错误：质控似乎未成功完成，请检查日志文件: ./logs/fastp.log"
    exit 1
fi
EOF

# 添加执行权限
chmod +x ./scripts/run_qc.sh

# 执行质控脚本
./scripts/run_qc.sh