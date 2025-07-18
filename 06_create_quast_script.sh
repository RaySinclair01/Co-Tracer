#!/bin/bash
# 文件名: create_quast_script.sh
# 用途: 创建质量评估脚本并执行

# 确保scripts目录存在
mkdir -p ./scripts

cat > ./scripts/run_quast.sh << 'EOF'
#!/bin/bash

# 定义输入文件
INPUT_FASTA="./input_data/contigs.fasta"

# 创建必要的目录
mkdir -p ./assembly
mkdir -p ./logs

# 检查输入文件
if [ ! -f "$INPUT_FASTA" ]; then
    echo "错误：输入文件 $INPUT_FASTA 不存在"
    echo "请确保组装步骤已成功完成并生成了contigs文件"
    exit 1
fi

# 删除已存在的QUAST报告目录(如果存在)
if [ -d "./assembly/quast_report" ]; then
    echo "发现已存在的QUAST报告目录，正在移除..."
    rm -rf ./assembly/quast_report
fi

# 使用QUAST评估组装质量
echo "评估组装质量..."
quast.py ${INPUT_FASTA} \
        -o ./assembly/quast_report \
        --threads 38 \
        > ./logs/quast.log 2>&1

# 检查QUAST是否成功运行
if [ -f "./assembly/quast_report/report.txt" ]; then
    echo "组装质量评估完成，结果如下："
    cat ./assembly/quast_report/report.txt
else
    echo "错误：QUAST评估似乎未成功完成"
    echo "请检查日志文件: ./logs/quast.log"
    exit 1
fi
EOF

# 添加执行权限
chmod +x ./scripts/run_quast.sh

# 执行质量评估脚本
./scripts/run_quast.sh