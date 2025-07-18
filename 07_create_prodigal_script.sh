#!/bin/bash
# 文件名: create_prodigal_script.sh
# 用途: 创建基因预测脚本并执行

# 确保scripts目录存在
mkdir -p ./scripts

cat > ./scripts/run_prodigal.sh << 'EOF'
#!/bin/bash

# 定义输入和输出文件路径
INPUT_FASTA="./input_data/contigs.fasta"
OUTPUT_DIR="./analysis_results"
OUTPUT_GFF="${OUTPUT_DIR}/predicted_genes.gff"
OUTPUT_FAA="${OUTPUT_DIR}/predicted_proteins.faa"
OUTPUT_GENES="${OUTPUT_DIR}/predicted_genes.fna"

# 确保输出目录和日志目录存在
mkdir -p ${OUTPUT_DIR}
mkdir -p ./logs

# 检查输入文件
if [ ! -f "$INPUT_FASTA" ]; then
    echo "错误：输入文件 $INPUT_FASTA 不存在"
    echo "请确保组装步骤已成功完成"
    exit 1
fi

# 执行基因预测
echo "开始基因预测..."
prodigal -i "${INPUT_FASTA}" \
         -o "${OUTPUT_GFF}" \
         -a "${OUTPUT_FAA}" \
         -d "${OUTPUT_GENES}" \
         -p meta -f gff -q \
         > ./logs/prodigal.log 2>&1

# 检查结果
if [ -f "$OUTPUT_FAA" ]; then
    echo "基因预测完成。结果文件:"
    ls -lh ${OUTPUT_DIR}

    # 简单统计
    echo "预测到的基因数量:"
    grep -c ">" ${OUTPUT_FAA}
else
    echo "错误：基因预测似乎未成功完成"
    echo "请检查日志文件: ./logs/prodigal.log"
    exit 1
fi
EOF

# 添加执行权限
chmod +x ./scripts/run_prodigal.sh

# 执行基因预测脚本
./scripts/run_prodigal.sh