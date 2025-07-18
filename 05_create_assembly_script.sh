#!/bin/bash
# 文件名: create_assembly_script.sh
# 用途: 创建组装脚本并执行

# 首先确保scripts目录存在
mkdir -p ./scripts

cat > ./scripts/run_assembly.sh << 'EOF'
#!/bin/bash

# 定义变量
QC_R1="./qc_data/L1HII0100002-TS_T1.R1.clean.fastq.gz"
QC_R2="./qc_data/L1HII0100002-TS_T1.R2.clean.fastq.gz"

# 创建必要的目录，但不创建megahit_out(让MEGAHIT自己创建)
mkdir -p ./assembly
mkdir -p ./logs
mkdir -p ./input_data

# 检查并移除已存在的输出目录
if [ -d "./assembly/megahit_out" ]; then
    echo "发现已存在的输出目录，正在移除..."
    rm -rf ./assembly/megahit_out
fi

# 双端数据组装 - 使用MEGAHIT
echo "开始双端数据组装..."
megahit -1 ${QC_R1} -2 ${QC_R2} \
        --num-cpu-threads 38 \
        --memory 0.9 \
        --min-contig-len 500 \
        --out-dir ./assembly/megahit_out \
        --out-prefix L1HII0100002-TS_T1 \
        > ./logs/megahit.log 2>&1

# 检查组装是否成功
if [ -f "./assembly/megahit_out/L1HII0100002-TS_T1.contigs.fa" ]; then
    echo "组装成功，找到contigs文件"
    ls -lh ./assembly/megahit_out/
    
    # 复制最终组装结果到输入目录
    cp ./assembly/megahit_out/L1HII0100002-TS_T1.contigs.fa ./input_data/contigs.fasta
    echo "组装完成，结果已复制到 ./input_data/contigs.fasta"
else
    echo "错误：组装失败，没有找到contigs文件"
    echo "查看日志文件了解详情: ./logs/megahit.log"
    exit 1
fi
EOF

# 添加执行权限
chmod +x ./scripts/run_assembly.sh

# 执行组装脚本
./scripts/run_assembly.sh