#!/bin/bash
# 文件名: create_arg_annotation_script.sh
# 用途: 创建ARG注释脚本并执行

cat > ./scripts/run_arg_annotation.sh << 'EOF'
#!/bin/bash

# 定义变量
OUTPUT_DIR="./analysis_results"
QUERY_FAA="${OUTPUT_DIR}/predicted_proteins.faa"
CARD_DB="$HOME/databases/card_protein_homolog_db"
OUTPUT_ARG_M8="${OUTPUT_DIR}/ARG_card_results.m8"

# 确保输出目录存在
mkdir -p "${OUTPUT_DIR}"
mkdir -p ./logs

# 检查数据库文件是否存在
echo "检查数据库文件..."
if [ -f "${CARD_DB}.dmnd" ]; then
    echo "找到数据库文件: ${CARD_DB}.dmnd"
else
    echo "错误: 未找到数据库文件: ${CARD_DB}.dmnd"
    echo "请确认数据库路径和文件名是否正确"
    exit 1
fi

# 检查输入文件
echo "检查输入蛋白文件..."
if [ -s "${QUERY_FAA}" ]; then
    echo "找到输入文件: ${QUERY_FAA}"
    echo "输入文件大小: $(ls -lh ${QUERY_FAA} | awk '{print $5}')"
    echo "序列数量: $(grep -c "^>" ${QUERY_FAA})"
else
    echo "错误: 未找到输入文件或文件为空: ${QUERY_FAA}"
    exit 1
fi

# 执行Diamond比对，使用更宽松的参数
echo "开始ARG注释(Diamond vs CARD)..."
diamond blastp \
    --db ${CARD_DB}.dmnd \
    --query ${QUERY_FAA} \
    --out ${OUTPUT_ARG_M8} \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp \
    --threads 38 \
    --evalue 1e-4 \
    --max-target-seqs 1 \
    --query-cover 70 \
    --id 70 \
    --sensitive \
    --block-size 6 \
    --index-chunks 2 \
    > ./logs/diamond.log 2>&1

# 检查结果
echo "Diamond完成。结果文件:"
ls -lh ${OUTPUT_ARG_M8}
wc -l ${OUTPUT_ARG_M8}

# 如果结果仍为空，提供额外信息
if [ ! -s "${OUTPUT_ARG_M8}" ]; then
    echo "警告: 结果文件为空，没有找到匹配项"
    echo "尝试检查日志文件获取更多信息:"
    tail -20 ./logs/diamond.log
    echo "可能需要进一步放宽搜索参数或检查数据库构建是否正确"
fi
EOF

# 添加执行权限
chmod +x ./scripts/run_arg_annotation.sh

# 执行ARG注释脚本
./scripts/run_arg_annotation.sh