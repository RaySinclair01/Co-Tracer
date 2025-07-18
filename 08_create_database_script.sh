#!/bin/bash
# 文件名: create_database_script.sh
# 用途: 创建数据库准备脚本并执行

cat > ./scripts/prepare_databases.sh << 'EOF'
#!/bin/bash

# 创建数据库目录
mkdir -p ~/databases
cd ~/databases

# 下载CARD数据库
echo "下载CARD数据库..."
wget https://card.mcmaster.ca/latest/data -O card_data.tar.bz2

# 解压CARD数据库
echo "解压CARD数据库..."
tar -xjf card_data.tar.bz2

# 构建Diamond数据库
echo "构建CARD Diamond数据库..."
diamond makedb --in protein_fasta_protein_homolog_model.fasta -d card_protein_homolog_db

# 下载Pfam-A HMM
echo "下载Pfam数据库..."
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/Pfam-A.hmm.gz -O Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

# 格式化Pfam数据库
echo "格式化Pfam数据库..."
hmmpress Pfam-A.hmm

echo "数据库准备完成"

# 返回原始工作目录
cd - > /dev/null
EOF

# 添加执行权限
chmod +x ./scripts/prepare_databases.sh

# 执行数据库准备脚本
./scripts/prepare_databases.sh