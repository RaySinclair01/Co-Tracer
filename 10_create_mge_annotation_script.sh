#!/bin/bash
# 文件名: create_mge_annotation_script.sh
# 用途: 创建MGE注释脚本并执行

cat > ./scripts/run_mge_annotation.sh << 'EOF'
#!/bin/bash

# 定义变量
OUTPUT_DIR="./analysis_results"
QUERY_FAA="${OUTPUT_DIR}/predicted_proteins.faa"
PFAM_HMM_DB="$HOME/databases/Pfam-A.hmm"
OUTPUT_PFAM_DOMTBL="${OUTPUT_DIR}/MGE_pfam_results.domtblout"

# 确保输出目录存在
mkdir -p "${OUTPUT_DIR}"
mkdir -p ./logs

# 检查数据库文件是否存在
echo "检查Pfam数据库文件..."
if [ -f "${PFAM_HMM_DB}" ]; then
    echo "找到数据库文件: ${PFAM_HMM_DB}"
else
    echo "错误: 未找到数据库文件: ${PFAM_HMM_DB}"
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

# 执行HMMER搜索
echo "开始MGE注释(HMMER vs Pfam)..."
hmmsearch \
    --domtblout ${OUTPUT_PFAM_DOMTBL} \
    --cut_ga \
    --cpu 16 \
    --noali \
    ${PFAM_HMM_DB} \
    ${QUERY_FAA} \
    > ./logs/hmmer_pfam.log 2>&1

# 检查结果
echo "HMMER Pfam搜索完成。结果文件:"
ls -lh ${OUTPUT_PFAM_DOMTBL}

# 计算匹配的条目数
MATCH_COUNT=$(grep -v "^#" ${OUTPUT_PFAM_DOMTBL} | wc -l)
echo "匹配条目数: ${MATCH_COUNT}"

# 如果结果为空或很少，提供额外信息
if [ "${MATCH_COUNT}" -lt 1 ]; then
    echo "警告: 结果文件没有匹配项"
    echo "尝试检查日志文件获取更多信息:"
    tail -20 ./logs/hmmer_pfam.log
    echo "可能需要检查Pfam数据库是否正确格式化(hmmpress)"
fi
EOF

# 添加执行权限
chmod +x ./scripts/run_mge_annotation.sh

# 执行MGE注释脚本
./scripts/run_mge_annotation.sh