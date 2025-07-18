#!/bin/bash
# 文件名: 11_process_annotations_and_colocalization.sh
# 用途: 创建数据处理脚本并执行，包含ARG/MGE处理、关联分析和共定位分析

cat > ./scripts/process_annotations_and_colocalization.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sqlite3
import os
import sys
import time
import logging
import re
import subprocess
from Bio import SeqIO
from collections import Counter

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("./logs/data_processing.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# 确保输出目录存在
os.makedirs("./logs", exist_ok=True)
os.makedirs("./analysis_results", exist_ok=True)

logger.info("开始数据处理与整合...")

# 检查文件是否存在
def check_file_exists(filepath, required=True):
    if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
        logger.info(f"文件存在且不为空: {filepath}")
        return True
    elif os.path.exists(filepath) and os.path.getsize(filepath) == 0:
        message = f"警告: 文件存在但为空: {filepath}"
        if required:
            logger.error(message)
            return False
        else:
            logger.warning(message)
            return False
    else:
        message = f"文件不存在: {filepath}"
        if required:
            logger.error(message)
            return False
        else:
            logger.warning(message)
            return False

# MGE相关的关键词列表 - 扩展版
MGE_KEYWORDS = [
    # 转座子相关
    'transpos', 'transposase', 'tnp', 'insertion', 'insertion sequence', 'is element',
    # 整合酶相关
    'integrase', 'integron', 'int', 'xis', 'excisionase', 
    # 噬菌体相关
    'phage', 'prophage', 'viral', 'virus', 'capsid', 'tail', 'head', 'baseplate',
    'portal', 'terminase', 'virion', 
    # 质粒相关
    'plasmid', 'replicon', 'replication', 'rep', 'replic', 'partition', 'par', 
    # 接合相关
    'conjugation', 'conjugative', 'conjugal', 'tra', 'trb', 'type iv', 
    # 分泌系统
    'secretion', 't4ss', 'type 4', 'type iv',
    # 移动相关
    'mobili', 'mobilization', 'mob', 'relaxase',
    # 重组相关
    'recomb', 'resolvase', 'invertase', 'flip', 'rci', 
    # 其他
    'transposon', 'retrotransposon', 'retron', 'reverse transcriptase', 'retrovir'
]

# 读取GFF文件，构建基因-坐标映射
def process_gff_file(filename):
    """读取GFF文件并构建ID映射"""
    logger.info(f"处理GFF文件: {filename}")
    
    if not check_file_exists(filename, required=True):
        logger.error("GFF文件不可用，终止处理")
        return pd.DataFrame(), {}, {}
    
    try:
        # 读取GFF文件
        gff_columns = ['contig', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
        gff_df = pd.read_csv(filename, sep='\t', comment='#', header=None, names=gff_columns)
        
        # 提取ID
        def extract_gff_id(attr_str):
            match = re.search(r'ID=([^;]+)', attr_str)
            if match:
                return match.group(1)
            return None
        
        gff_df['gene_id'] = gff_df['attributes'].apply(extract_gff_id)
        
        # 构建k141_style_id
        def build_k141_id(row):
            if not pd.isna(row['gene_id']) and '_' in str(row['gene_id']):
                contig_base = row['contig']  # 例如 k141_45964
                gene_num = row['gene_id'].split('_')[1]  # 例如 从1_2中获取2
                return f"{contig_base}_{gene_num}"
            return None
        
        gff_df['k141_style_id'] = gff_df.apply(build_k141_id, axis=1)
        
        # 创建ID映射
        k141_to_coords = {}
        gff_id_to_coords = {}
        
        for _, row in gff_df.iterrows():
            coords = {
                'contig': row['contig'],
                'start': row['start'],
                'end': row['end'],
                'strand': row['strand'],
                'gene_id': row['gene_id']
            }
            
            if not pd.isna(row['k141_style_id']):
                k141_to_coords[row['k141_style_id']] = coords
            
            if not pd.isna(row['gene_id']):
                gff_id_to_coords[row['gene_id']] = coords
        
        logger.info(f"GFF文件中的条目数: {len(gff_df)}")
        logger.info(f"构建了 {len(k141_to_coords)} 个k141风格ID映射")
        logger.info(f"构建了 {len(gff_id_to_coords)} 个GFF ID映射")
        
        # 显示一些示例
        logger.info("GFF文件示例ID:")
        for idx, row in gff_df.head(3).iterrows():
            logger.info(f"原始ID: {row['gene_id']}, 转换后ID: {row['k141_style_id']}, Contig: {row['contig']}")
        
        return gff_df, k141_to_coords, gff_id_to_coords
    
    except Exception as e:
        logger.error(f"处理GFF文件时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return pd.DataFrame(), {}, {}

# 处理ARG数据
def process_arg_data(arg_file, k141_to_coords, gff_id_to_coords):
    """处理ARG注释数据并关联坐标"""
    logger.info(f"处理ARG注释结果: {arg_file}")
    
    if not check_file_exists(arg_file, required=False):
        logger.warning("ARG文件不存在或为空，跳过ARG处理")
        return pd.DataFrame()
    
    try:
        arg_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                    'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                    'evalue', 'bitscore', 'qcovhsp']
                    
        arg_df = pd.read_csv(arg_file, sep='\t', names=arg_cols)
        
        # 显示ARG查询ID样本
        if not arg_df.empty:
            arg_id_samples = arg_df['qseqid'].head(5).tolist()
            logger.info(f"ARG查询ID样本: {', '.join(map(str, arg_id_samples))}")
        
        # 检查是否有数据
        if len(arg_df) == 0:
            logger.warning("ARG文件为空，跳过ARG处理")
            return pd.DataFrame()
        
        # 过滤ARG结果
        arg_filtered = arg_df[(arg_df['evalue'] <= 1e-10) & 
                           (arg_df['bitscore'] >= 60) & 
                           (arg_df['pident'] >= 80) & 
                           (arg_df['qcovhsp'] >= 80)].copy()

        logger.info(f"过滤后的ARG注释: {len(arg_filtered)} / {len(arg_df)}")
        
        # 关联坐标
        arg_with_coords = pd.DataFrame()
        match_count = 0
        
        # 1. 通过k141_style_id匹配
        for idx, row in arg_filtered.iterrows():
            if row['qseqid'] in k141_to_coords:
                coords = k141_to_coords[row['qseqid']]
                arg_filtered.loc[idx, 'contig'] = coords['contig']
                arg_filtered.loc[idx, 'start'] = coords['start']
                arg_filtered.loc[idx, 'end'] = coords['end']
                arg_filtered.loc[idx, 'strand'] = coords['strand']
                arg_filtered.loc[idx, 'gene_id'] = coords['gene_id']
                match_count += 1
        
        logger.info(f"通过k141 ID匹配的ARG条目: {match_count} / {len(arg_filtered)}")
        
        # 清理无法关联坐标的条目
        arg_with_coords = arg_filtered.dropna(subset=['contig', 'start', 'end'])
        logger.info(f"最终成功关联坐标的ARG: {len(arg_with_coords)} / {len(arg_filtered)}")
        
        # 如果有匹配的ARG，处理CARD信息
        if not arg_with_coords.empty:
            # 从CARD ID中提取ARG信息
            def parse_card_id(sseqid):
                parts = str(sseqid).split('|')
                
                # 默认值
                aro = 'ARO:unknown'
                gene_name = 'unknown'
                
                # 查找ARO部分
                for i, part in enumerate(parts):
                    if 'ARO:' in part:
                        aro = part
                        # 如果ARO后面还有部分，那就是亚类
                        if i + 1 < len(parts):
                            gene_name = parts[i + 1]  # 提取ARO后面的部分作为亚类
                        break
                
                return {'aro': aro, 'gene_name': gene_name}

            # 应用ARG解析
            arg_info = arg_with_coords['sseqid'].apply(parse_card_id)
            arg_with_coords.loc[:, 'aro'] = arg_info.apply(lambda x: x['aro'])
            arg_with_coords.loc[:, 'gene_name'] = arg_info.apply(lambda x: x['gene_name'])
            arg_with_coords.loc[:, 'arg_family'] = arg_with_coords['gene_name']
            
            # 添加ARG亚类信息
            arg_with_coords.loc[:, 'arg_subtype'] = arg_with_coords['gene_name']
            
            # 确保有qseqid作为基因ID
            if 'gene_id' not in arg_with_coords.columns or arg_with_coords['gene_id'].isna().any():
                arg_with_coords.loc[:, 'gene_id'] = arg_with_coords['qseqid']
        
        return arg_with_coords
    
    except Exception as e:
        logger.error(f"处理ARG数据时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return pd.DataFrame()

# 创建ARG数据库
def create_arg_database(arg_data, output_file):
    """创建ARG数据库"""
    logger.info(f"创建ARG SQLite数据库: {output_file}")
    
    if arg_data.empty:
        logger.warning("没有ARG数据可写入数据库")
        return False
    
    try:
        # 如果数据库已存在，先删除
        if os.path.exists(output_file):
            os.remove(output_file)
        
        # 创建SQLite连接
        conn = sqlite3.connect(output_file)
        cursor = conn.cursor()
        
        # 创建ARG表
        cursor.execute('''
        CREATE TABLE arg_hits (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_id TEXT,
            contig TEXT,
            start INTEGER,
            end INTEGER,
            strand TEXT,
            aro TEXT,
            gene_name TEXT,
            arg_family TEXT,
            arg_subtype TEXT,
            pident REAL,
            bitscore REAL,
            evalue REAL,
            qcovhsp REAL
        )
        ''')
        
        # 创建索引
        cursor.execute('CREATE INDEX idx_arg_contig ON arg_hits(contig)')
        cursor.execute('CREATE INDEX idx_arg_family ON arg_hits(arg_family)')
        cursor.execute('CREATE INDEX idx_arg_subtype ON arg_hits(arg_subtype)')
        
        # 准备ARG数据
        arg_columns = ['gene_id', 'contig', 'start', 'end', 'strand',
                      'aro', 'gene_name', 'arg_family', 'arg_subtype', 'pident', 
                      'bitscore', 'evalue', 'qcovhsp']
        
        # 确保所有必需列都存在
        for col in arg_columns:
            if col not in arg_data.columns:
                logger.warning(f"列 {col} 不存在于ARG数据中，使用NaN填充")
                arg_data[col] = np.nan
        
        # 插入数据
        insert_data = []
        for _, row in arg_data.iterrows():
            insert_data.append((
                row['gene_id'], row['contig'], row['start'], row['end'],
                row['strand'], row['aro'], row['gene_name'], row['arg_family'],
                row['arg_subtype'], row['pident'], row['bitscore'], row['evalue'], row['qcovhsp']
            ))
        
        cursor.executemany('''
        INSERT INTO arg_hits 
        (gene_id, contig, start, end, strand, aro, gene_name, arg_family, arg_subtype, pident, bitscore, evalue, qcovhsp)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', insert_data)
        
        # 提交并关闭
        conn.commit()
        
        # 验证插入
        cursor.execute("SELECT COUNT(*) FROM arg_hits")
        count = cursor.fetchone()[0]
        logger.info(f"已成功写入 {count} 条ARG记录到数据库")
        
        conn.close()
        return True
    
    except Exception as e:
        logger.error(f"创建ARG数据库时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

# 解析HMMER domtblout文件
def parse_domtblout(filename, batch_size=100000):
    """分批解析HMMER domtblout文件，返回生成器"""
    logger.info(f"开始解析domtblout文件: {filename}")
    
    if not check_file_exists(filename, required=False):
        logger.warning("MGE文件不存在或为空，跳过MGE处理")
        return
    
    try:
        records = []
        line_count = 0
        valid_count = 0
        batch_num = 1
        pfam_counts = Counter()
        pfam_names = {}
        mge_related_pfams = set()
        
        with open(filename, 'r') as f:
            for line in f:
                line_count += 1
                
                if line.startswith('#') or not line.strip():
                    continue
                
                try:
                    # 按空格分割，限制最大分割数为22
                    parts = line.strip().split(None, 22)
                    
                    if len(parts) < 12:  # 确保至少有必需的列
                        continue
                    
                    # 提取基本信息
                    target_name = parts[0]      # 序列ID
                    hmm_name = parts[3]         # Pfam域名称
                    accession = parts[4]        # Pfam ID
                    e_value = float(parts[6])   # 全序列E值
                    score = float(parts[7])     # 位得分
                    i_evalue = float(parts[11]) # 独立E值
                    
                    # 提取描述中的GFF ID
                    gff_id = None
                    if len(parts) > 22:
                        id_match = re.search(r'ID=([^;]+)', parts[22])
                        if id_match:
                            gff_id = id_match.group(1)
                    
                    # 收集Pfam统计
                    pfam_counts[accession] += 1
                    pfam_names[accession] = hmm_name
                    
                    # 检查是否是MGE相关域
                    hmm_name_lower = hmm_name.lower()
                    is_mge_related = any(keyword in hmm_name_lower for keyword in MGE_KEYWORDS)
                    
                    if is_mge_related:
                        mge_related_pfams.add(accession)
                    
                    # 分类MGE类型
                    mge_category = 'Other MGE'
                    if any(kw in hmm_name_lower for kw in ['integrase', 'integron', 'int']):
                        mge_category = 'Integrase'
                    elif any(kw in hmm_name_lower for kw in ['transpos', 'tnp', 'insertion']):
                        mge_category = 'Transposase'
                    elif any(kw in hmm_name_lower for kw in ['phage', 'virus', 'viral', 'capsid', 'tail']):
                        mge_category = 'Phage'
                    elif any(kw in hmm_name_lower for kw in ['plasmid', 'replicon', 'rep']):
                        mge_category = 'Plasmid'
                    elif any(kw in hmm_name_lower for kw in ['conjug', 'tra', 'trb', 'pili']):
                        mge_category = 'Conjugation'
                    elif any(kw in hmm_name_lower for kw in ['mobili', 'mob', 'relaxase']):
                        mge_category = 'Mobilization'
                    elif any(kw in hmm_name_lower for kw in ['recomb', 'resolvase', 'invertase']):
                        mge_category = 'Recombination'
                    
                    # 添加亚类信息
                    mge_subtype = hmm_name
                    
                    # 高置信度匹配：高分 + 低E值 或 MGE相关
                    high_confidence = score >= 25 and i_evalue <= 1e-5
                    
                    # 如果是高置信度匹配或MGE相关
                    if high_confidence or is_mge_related:
                        records.append({
                            'target_name': target_name,
                            'accession': accession,
                            'hmm_name': hmm_name,
                            'e_value': e_value,
                            'score': score,
                            'i_evalue': i_evalue,
                            'gff_id': gff_id,
                            'mge_category': mge_category,
                            'mge_family': f"{accession}:{hmm_name}",
                            'mge_subtype': mge_subtype
                        })
                        
                        valid_count += 1
                        
                        # 返回一个批次的数据
                        if len(records) >= batch_size:
                            if batch_num == 1:
                                # 首批次时，同时返回Pfam统计信息
                                logger.info(f"已处理 {line_count} 行，返回第 {batch_num} 批记录（{len(records)} 条）")
                                common_pfams = pfam_counts.most_common(20)
                                logger.info("最常见的20个Pfam域:")
                                for pfam, count in common_pfams:
                                    is_mge = "（可能是MGE）" if pfam in mge_related_pfams else ""
                                    logger.info(f"  {pfam} ({pfam_names[pfam]}): {count}次 {is_mge}")
                                
                                logger.info(f"基于关键词识别了 {len(mge_related_pfams)} 个MGE相关域")
                            else:
                                logger.info(f"已处理 {line_count} 行，返回第 {batch_num} 批记录（{len(records)} 条）")
                            
                            yield records
                            records = []
                            batch_num += 1
                
                except (ValueError, IndexError) as e:
                    if line_count % 100000 == 0:
                        logger.warning(f"解析第 {line_count} 行时出错: {str(e)}")
                    continue
                
                # 每处理100万行输出一次日志
                if line_count % 1000000 == 0:
                    logger.info(f"已处理 {line_count} 行，找到 {valid_count} 条有效记录")
        
        # 返回最后一批数据
        if records:
            logger.info(f"已处理完全部 {line_count} 行，返回最后一批记录（{len(records)} 条）")
            yield records
        
        logger.info(f"共解析 {line_count} 行，找到 {valid_count} 条有效记录")
        logger.info(f"基于关键词识别了共 {len(mge_related_pfams)} 个MGE相关域")
    
    except Exception as e:
        logger.error(f"解析domtblout文件时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        yield []

# 创建MGE数据库
def create_mge_database(mge_file, k141_to_coords, gff_id_to_coords, output_file):
    """分批处理并创建MGE数据库"""
    logger.info(f"创建MGE数据库: {output_file}")
    
    try:
        # 如果数据库已存在，先删除
        if os.path.exists(output_file):
            logger.info(f"删除现有数据库: {output_file}")
            os.remove(output_file)
        
        # 创建数据库连接
        conn = sqlite3.connect(output_file)
        cursor = conn.cursor()
        
        # 创建MGE表
        cursor.execute('''
        CREATE TABLE mge_hits (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_id TEXT,
            target_name TEXT,
            contig TEXT,
            start INTEGER,
            end INTEGER,
            strand TEXT,
            pfam_accession TEXT,
            pfam_name TEXT,
            mge_category TEXT,
            mge_family TEXT,
            mge_subtype TEXT,
            score REAL,
            i_evalue REAL
        )
        ''')
        
        # 创建索引
        cursor.execute('CREATE INDEX idx_mge_contig ON mge_hits(contig)')
        cursor.execute('CREATE INDEX idx_mge_category ON mge_hits(mge_category)')
        cursor.execute('CREATE INDEX idx_mge_subtype ON mge_hits(mge_subtype)')
        
        # 提交初始架构
        conn.commit()
        
        # 分批处理MGE数据
        total_inserted = 0
        batch_num = 1
        match_methods = Counter()
        
        for batch in parse_domtblout(mge_file, batch_size=50000):
            logger.info(f"处理第 {batch_num} 批次MGE数据（{len(batch)} 条记录）")
            
            # 关联坐标
            records_with_coords = []
            
            for record in batch:
                # 方法1：通过k141_style_id匹配
                if record['target_name'] in k141_to_coords:
                    coords = k141_to_coords[record['target_name']]
                    record.update(coords)
                    records_with_coords.append(record)
                    match_methods['k141_id'] += 1
                    continue
                
                # 方法2：通过gff_id匹配
                if record['gff_id'] and record['gff_id'] in gff_id_to_coords:
                    coords = gff_id_to_coords[record['gff_id']]
                    record.update(coords)
                    records_with_coords.append(record)
                    match_methods['gff_id'] += 1
                    continue
            
            logger.info(f"批次 {batch_num}: 成功关联 {len(records_with_coords)} / {len(batch)} 条记录")
            
            # 插入数据
            if records_with_coords:
                insert_data = []
                for record in records_with_coords:
                    insert_data.append((
                        record.get('gene_id', ''),
                        record.get('target_name', ''),
                        record.get('contig', ''),
                        record.get('start', 0),
                        record.get('end', 0),
                        record.get('strand', ''),
                        record.get('accession', ''),
                        record.get('hmm_name', ''),
                        record.get('mge_category', ''),
                        record.get('mge_family', ''),
                        record.get('mge_subtype', ''),
                        record.get('score', 0.0),
                        record.get('i_evalue', 0.0)
                    ))
                
                cursor.executemany('''
                INSERT INTO mge_hits 
                (gene_id, target_name, contig, start, end, strand,
                 pfam_accession, pfam_name, mge_category, mge_family, mge_subtype, score, i_evalue)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', insert_data)
                
                # 提交本批次
                conn.commit()
                
                # 验证插入
                cursor.execute("SELECT COUNT(*) FROM mge_hits")
                count = cursor.fetchone()[0]
                
                logger.info(f"批次 {batch_num}: 插入 {len(insert_data)} 条记录，数据库中现有 {count} 条记录")
                total_inserted += len(insert_data)
            
            batch_num += 1
        
        # 最终验证
        cursor.execute("SELECT COUNT(*) FROM mge_hits")
        final_count = cursor.fetchone()[0]
        
        cursor.execute("SELECT mge_category, COUNT(*) FROM mge_hits GROUP BY mge_category")
        category_counts = cursor.fetchall()
        
        logger.info(f"MGE数据库创建完成，共插入 {total_inserted} 条记录，验证有 {final_count} 条记录")
        logger.info("MGE分类统计:")
        for category, count in category_counts:
            logger.info(f"  {category}: {count}")
        
        logger.info("坐标匹配方法统计:")
        for method, count in match_methods.most_common():
            logger.info(f"  {method}: {count}")
        
        conn.close()
        
        # 确认数据库创建成功
        if final_count > 0:
            logger.info("MGE数据库创建成功")
            return True
        else:
            logger.error("MGE数据库创建失败，没有记录被插入")
            return False
    
    except Exception as e:
        logger.error(f"创建MGE数据库时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

# 执行contig分类学注释
def perform_taxonomy_analysis(common_contigs):
    """对包含ARG和MGE的contig进行分类学注释"""
    logger.info("开始对共定位contig进行分类学注释...")
    
    # 检查是否有共定位的contig
    if not common_contigs:
        logger.warning("没有共定位的contig，跳过分类学分析")
        return {}
    
    # 定义文件路径
    input_contigs = "./input_data-HDS/contigs.fasta"  # 修正的路径
    colocalized_contigs_list = "./analysis_results-HDS/colocalized_contigs.list"
    colocalized_contigs_fasta = "./analysis_results-HDS/colocalized_contigs.fasta"
    output_prefix = "./analysis_results-HDS/taxonomy"
    tax_csv = "./analysis_results-HDS/contig_taxonomy.csv"
    
    # 检查输入contigs文件
    if not check_file_exists(input_contigs, required=False):
        logger.error(f"输入contigs文件不存在或为空: {input_contigs}")
        return {}
    
    try:
        # 保存共定位contig列表
        with open(colocalized_contigs_list, 'w') as f:
            for contig in common_contigs:
                f.write(f"{contig}\n")
        
        logger.info(f"已保存 {len(common_contigs)} 个共定位contig ID")
        
        # 提取共定位contig序列
        extracted_count = 0
        try:
            with open(colocalized_contigs_fasta, 'w') as out_f:
                for record in SeqIO.parse(input_contigs, "fasta"):
                    if record.id in common_contigs:
                        SeqIO.write(record, out_f, "fasta")
                        extracted_count += 1
            
            logger.info(f"已提取 {extracted_count} 个contig序列到文件: {colocalized_contigs_fasta}")
        except Exception as e:
            logger.error(f"提取contig序列时出错: {str(e)}")
            return {}
        
        # 检查提取的序列文件
        if not check_file_exists(colocalized_contigs_fasta, required=True):
            logger.error("提取的contig序列文件为空或不存在")
            return {}
        
        # 检查CAT/BAT是否可用
        cat_path = os.path.expanduser("/mnt/d/CAT_database")
        cat_tax_path = os.path.expanduser("/mnt/d/CAT_taxonomy")
        
        if not os.path.exists(cat_path) or not os.path.exists(cat_tax_path):
            logger.warning("CAT/BAT数据库不存在，尝试寻找其他方式")
            # 这里可以添加其他分类器的支持
            return {}
        
        # 运行CAT进行分类学注释
        try:
            logger.info("运行CAT进行分类学注释...")
            cat_cmd = [
                "CAT_pack", "contigs",
                "-c", colocalized_contigs_fasta,
                "-d", cat_path,
                "-t", cat_tax_path,
                "-o", output_prefix,
                "--force",
                "-n", "38",
                "--sensitive"
            ]
            
            logger.info(f"执行命令: {' '.join(cat_cmd)}")
            cat_process = subprocess.run(cat_cmd, capture_output=True, text=True)
            
            if cat_process.returncode != 0:
                logger.error(f"CAT执行失败: {cat_process.stderr}")
                return {}
            
            # 添加分类学名称
            names_cmd = [
                "CAT_pack", "add_names",
                "-i", f"{output_prefix}.contig2classification.txt",
                "-o", f"{output_prefix}.contig2classification.names.txt",
                "-t", cat_tax_path
            ]
            
            logger.info(f"执行命令: {' '.join(names_cmd)}")
            names_process = subprocess.run(names_cmd, capture_output=True, text=True)
            
            if names_process.returncode != 0:
                logger.error(f"CAT add_names执行失败: {names_process.stderr}")
                return {}
            
            # 解析分类学结果
            tax_file = f"{output_prefix}.contig2classification.names.txt"
            
            if not check_file_exists(tax_file, required=True):
                logger.error(f"分类学结果文件不存在或为空: {tax_file}")
                return {}
            
            tax_data = []
            with open(tax_file, 'r') as f:
                next(f)  # 跳过标题行
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 6:  # 确保有足够的列
                        continue
                    
                    contig_id = parts[0]
                    classification = parts[1]
                    
                    # 提取门和属信息
                    phylum = "unknown"
                    genus = "unknown"
                    
                    # 完整谱系信息从第6列开始
                    full_lineage_parts = parts[5:]
                    for item in full_lineage_parts:
                        if "(phylum):" in item:
                            phylum = item.split('(phylum):')[0].strip()
                        if "(genus):" in item:
                            genus = item.split('(genus):')[0].strip()
                    
                    tax_data.append({
                        'contig': contig_id,
                        'classification': classification,
                        'phylum': phylum,
                        'genus': genus,
                        'full_lineage': ' '.join(full_lineage_parts)
                    })
            
            # 创建DataFrame并保存
            if tax_data:
                tax_df = pd.DataFrame(tax_data)
                tax_df.to_csv(tax_csv, index=False)
                logger.info(f"已创建分类学CSV文件，包含 {len(tax_df)} 条contig的分类信息")
                
                # 创建分类学字典返回
                tax_dict = dict(zip(tax_df['contig'], zip(tax_df['phylum'], tax_df['genus'])))
                return tax_dict
            else:
                logger.warning("没有分类学数据可用")
                return {}
            
        except Exception as e:
            logger.error(f"运行CAT分类学注释时出错: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            return {}
            
    except Exception as e:
        logger.error(f"执行分类学分析时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return {}

# 创建ARG与MGE关联分析
def create_arg_mge_association(arg_db, mge_db, output_db):
    """创建ARG与MGE距离关联分析"""
    logger.info("创建ARG与MGE距离关联分析...")
    
    try:
        # 检查数据库文件是否存在
        if not os.path.exists(arg_db):
            logger.error(f"ARG数据库文件不存在: {arg_db}")
            return False, [], set()
        
        if not os.path.exists(mge_db):
            logger.error(f"MGE数据库文件不存在: {mge_db}")
            return False, [], set()
        
        # 使用SQLite直接连接
        arg_conn = sqlite3.connect(arg_db)
        mge_conn = sqlite3.connect(mge_db)
        
        # 读取ARG数据
        arg_df = pd.read_sql("SELECT * FROM arg_hits", arg_conn)
        logger.info(f"从ARG数据库读取了 {len(arg_df)} 条记录")
        
        # 读取MGE数据
        mge_df = pd.read_sql("SELECT * FROM mge_hits", mge_conn)
        logger.info(f"从MGE数据库读取了 {len(mge_df)} 条记录")
        
        # 关闭连接
        arg_conn.close()
        mge_conn.close()
        
        # 确保有数据可用
        if arg_df.empty:
            logger.error("ARG数据为空，无法创建关联")
            return False, [], set()
        
        if mge_df.empty:
            logger.error("MGE数据为空，无法创建关联")
            return False, [], set()
        
        # 计算ARG与MGE的距离
        logger.info("计算ARG与MGE的距离关联...")
        arg_mge_pairs = []
        
        # 收集共定位的contig
        common_contigs = set()
        
        for _, arg in arg_df.iterrows():
            arg_contig = arg['contig']
            arg_start = int(arg['start'])
            arg_end = int(arg['end'])
            
            # 查找同一contig上的MGE
            same_contig_mges = mge_df[mge_df['contig'] == arg_contig]
            
            if not same_contig_mges.empty:
                common_contigs.add(arg_contig)
            
            for _, mge in same_contig_mges.iterrows():
                mge_start = int(mge['start'])
                mge_end = int(mge['end'])
                
                # 计算距离（如果重叠，距离为0）
                if arg_end < mge_start:  # ARG在MGE之前
                    distance = mge_start - arg_end
                elif mge_end < arg_start:  # MGE在ARG之前
                    distance = arg_start - mge_end
                else:  # 重叠
                    distance = 0
                
                # 添加到对列表
                arg_mge_pairs.append({
                    'arg_id': arg['gene_id'],
                    'arg_family': arg.get('arg_family', ''),
                    'arg_subtype': arg.get('arg_subtype', ''),
                    'mge_id': mge['gene_id'],
                    'mge_category': mge['mge_category'],
                    'mge_family': mge['mge_family'],
                    'mge_subtype': mge.get('mge_subtype', ''),
                    'contig': arg_contig,
                    'distance': distance,
                    'arg_start': arg_start,
                    'arg_end': arg_end,
                    'mge_start': mge_start,
                    'mge_end': mge_end
                })
        
        # 创建关联分析数据库
        if arg_mge_pairs:
            logger.info(f"共找到 {len(arg_mge_pairs)} 个ARG-MGE关联，涉及 {len(common_contigs)} 个共定位contig")
            
            # 如果数据库已存在，先删除
            if os.path.exists(output_db):
                os.remove(output_db)
            
            # 创建新数据库
            assoc_conn = sqlite3.connect(output_db)
            cursor = assoc_conn.cursor()
            
            # 创建关联表
            cursor.execute('''
            CREATE TABLE arg_mge_associations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                arg_id TEXT,
                arg_family TEXT,
                arg_subtype TEXT,
                mge_id TEXT,
                mge_category TEXT,
                mge_family TEXT,
                mge_subtype TEXT,
                contig TEXT,
                distance INTEGER,
                arg_start INTEGER,
                arg_end INTEGER,
                mge_start INTEGER,
                mge_end INTEGER,
                phylum TEXT,
                genus TEXT
            )
            ''')
            
            # 插入数据
            insert_data = []
            for pair in arg_mge_pairs:
                insert_data.append((
                    pair['arg_id'], pair['arg_family'], pair['arg_subtype'], 
                    pair['mge_id'], pair['mge_category'], pair['mge_family'], pair['mge_subtype'],
                    pair['contig'], pair['distance'], 
                    pair['arg_start'], pair['arg_end'], pair['mge_start'], pair['mge_end'],
                    'unknown', 'unknown'  # 先用unknown占位，后面会更新
                ))
            
            # 分批插入，避免一次插入太多数据
            batch_size = 10000
            for i in range(0, len(insert_data), batch_size):
                batch = insert_data[i:i + batch_size]
                cursor.executemany('''
                INSERT INTO arg_mge_associations 
                (arg_id, arg_family, arg_subtype, mge_id, mge_category, mge_family, mge_subtype, 
                 contig, distance, arg_start, arg_end, mge_start, mge_end, phylum, genus)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', batch)
                assoc_conn.commit()
            
            # 创建统计表
            cursor.execute('''
            CREATE TABLE distance_statistics (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                arg_family TEXT,
                mge_category TEXT,
                min REAL,
                max REAL,
                mean REAL,
                count INTEGER
            )
            ''')
            
            # 将数据转换为DataFrame进行统计
            pairs_df = pd.DataFrame(arg_mge_pairs)
            distance_stats = pairs_df.groupby(['arg_family', 'mge_category'])['distance'].agg(['min', 'max', 'mean', 'count']).reset_index()
            
            # 插入统计数据
            for _, row in distance_stats.iterrows():
                cursor.execute('''
                INSERT INTO distance_statistics 
                (arg_family, mge_category, min, max, mean, count)
                VALUES (?, ?, ?, ?, ?, ?)
                ''', (row['arg_family'], row['mge_category'], row['min'], row['max'], row['mean'], row['count']))
            
            assoc_conn.commit()
            
            # 验证数据
            cursor.execute("SELECT COUNT(*) FROM arg_mge_associations")
            assoc_count = cursor.fetchone()[0]
            
            cursor.execute("SELECT COUNT(*) FROM distance_statistics")
            stats_count = cursor.fetchone()[0]
            
            logger.info(f"关联数据库创建完成:")
            logger.info(f"  - arg_mge_associations表: {assoc_count}条记录")
            logger.info(f"  - distance_statistics表: {stats_count}条记录")
            
            # 统计信息
            overlap_count = len(pairs_df[pairs_df['distance'] == 0])
            close_count = len(pairs_df[pairs_df['distance'] <= 10000])
            
            logger.info(f"重叠的ARG-MGE对: {overlap_count} ({overlap_count/len(pairs_df)*100:.1f}%)")
            logger.info(f"距离≤10kb的ARG-MGE对: {close_count} ({close_count/len(pairs_df)*100:.1f}%)")
            
            assoc_conn.close()
            return True, pairs_df, common_contigs
        else:
            logger.warning("没有找到ARG-MGE关联，跳过创建关联数据库")
            return False, pd.DataFrame(), set()
    
    except Exception as e:
        logger.error(f"创建ARG-MGE关联时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False, pd.DataFrame(), set()

# 更新分类学信息到关联数据库
def update_taxonomy_in_database(tax_dict, output_db):
    """更新关联数据库中的分类学信息"""
    if not tax_dict:
        logger.warning("没有分类学信息可用，跳过更新")
        return False
    
    if not os.path.exists(output_db):
        logger.error(f"关联数据库不存在: {output_db}")
        return False
    
    try:
        # 连接数据库
        conn = sqlite3.connect(output_db)
        cursor = conn.cursor()
        
        # 为每个contig更新分类学信息
        update_count = 0
        for contig, (phylum, genus) in tax_dict.items():
            cursor.execute('''
            UPDATE arg_mge_associations
            SET phylum = ?, genus = ?
            WHERE contig = ?
            ''', (phylum, genus, contig))
            
            update_count += cursor.rowcount
        
        # 提交更改
        conn.commit()
        
        # 验证更新
        cursor.execute("SELECT COUNT(*) FROM arg_mge_associations WHERE phylum != 'unknown'")
        classified_count = cursor.fetchone()[0]
        
        logger.info(f"更新了 {update_count} 条记录的分类学信息，现有 {classified_count} 条记录有分类信息")
        
        # 关闭连接
        conn.close()
        return True
    
    except Exception as e:
        logger.error(f"更新分类学信息时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

# 生成共定位分析CSV报告
def generate_colocalization_reports(pairs_df, tax_dict):
    """生成共定位分析的CSV报告，包含亚类和分类学信息"""
    logger.info("生成共定位分析CSV报告...")
    
    if pairs_df.empty:
        logger.warning("关联对数据为空，创建空的共定位报告")
        pd.DataFrame().to_csv('./analysis_results-HDS/co_localization_details.csv', index=False)
        pd.DataFrame(columns=['arg_family', 'mge_category', 'count']).to_csv('./analysis_results-HDS/co_localization_summary_counts.csv', index=False)
        pd.DataFrame(columns=['arg_family', 'mge_family', 'count']).to_csv('./analysis_results-HDS/co_localization_detailed_counts.csv', index=False)
        pd.DataFrame(columns=['arg_subtype', 'mge_subtype', 'count']).to_csv('./analysis_results-HDS/co_localization_subtype_counts.csv', index=False)
        logger.info("已创建空的共定位报告文件")
        return
    
    try:
        # 设置共定位阈值
        DISTANCE_THRESHOLD = 10000  # 10kb
        logger.info(f"使用共定位距离阈值: {DISTANCE_THRESHOLD}bp")
        
        # 筛选共定位的对
        coloc_df = pairs_df[pairs_df['distance'] <= DISTANCE_THRESHOLD].copy()
        logger.info(f"发现 {len(coloc_df)} 条共定位事件（距离 ≤ {DISTANCE_THRESHOLD}bp）")
        
        # 添加分类学信息
        if tax_dict:
            logger.info("添加分类学信息到共定位结果...")
            
            # 添加门和属列
            coloc_df['phylum'] = coloc_df['contig'].apply(lambda x: tax_dict.get(x, ('unknown', 'unknown'))[0])
            coloc_df['genus'] = coloc_df['contig'].apply(lambda x: tax_dict.get(x, ('unknown', 'unknown'))[1])
            
            # 统计分类学信息
            phylum_count = coloc_df['phylum'].nunique()
            genus_count = coloc_df['genus'].nunique()
            classified_count = len(coloc_df[coloc_df['phylum'] != 'unknown'])
            
            logger.info(f"分类学统计: {phylum_count}个门, {genus_count}个属")
            logger.info(f"有分类信息的共定位事件: {classified_count}/{len(coloc_df)} ({classified_count/len(coloc_df)*100:.1f}%)")
        else:
            logger.warning("没有分类学信息可用，使用'unknown'作为默认值")
            coloc_df['phylum'] = 'unknown'
            coloc_df['genus'] = 'unknown'
        
        if len(coloc_df) == 0:
            logger.warning("没有符合条件的共定位事件")
            # 创建空结果文件
            coloc_df.to_csv('./analysis_results-HDS/co_localization_details.csv', index=False)
            pd.DataFrame(columns=['arg_family', 'mge_category', 'count']).to_csv('./analysis_results-HDS/co_localization_summary_counts.csv', index=False)
            pd.DataFrame(columns=['arg_family', 'mge_family', 'count']).to_csv('./analysis_results-HDS/co_localization_detailed_counts.csv', index=False)
            pd.DataFrame(columns=['arg_subtype', 'mge_subtype', 'count']).to_csv('./analysis_results-HDS/co_localization_subtype_counts.csv', index=False)
            logger.info("已创建空的共定位报告文件")
            return
        
        # 创建统计信息
        # 按ARG家族和MGE类别统计
        category_counts = coloc_df.groupby(['arg_family', 'mge_category']).size().reset_index(name='count')
        
        # 按ARG家族和MGE具体家族统计
        family_counts = coloc_df.groupby(['arg_family', 'mge_family']).size().reset_index(name='count')
        
        # 按ARG亚类和MGE亚类统计
        subtype_counts = coloc_df.groupby(['arg_subtype', 'mge_subtype']).size().reset_index(name='count')
        
        # 保存结果
        coloc_df.to_csv('./analysis_results-HDS/co_localization_details.csv', index=False)
        category_counts.to_csv('./analysis_results-HDS/co_localization_summary_counts.csv', index=False)
        family_counts.to_csv('./analysis_results-HDS/co_localization_detailed_counts.csv', index=False)
        subtype_counts.to_csv('./analysis_results-HDS/co_localization_subtype_counts.csv', index=False)
        
        # 提供统计摘要
        overlap_count = len(coloc_df[coloc_df['distance'] == 0])
        logger.info(f"共定位统计摘要:")
        logger.info(f"  - 总共定位事件: {len(coloc_df)}条")
        logger.info(f"  - 重叠ARG-MGE对: {overlap_count}条 ({overlap_count/len(coloc_df)*100:.1f}%)")
        
        # 距离分布
        logger.info("  - 共定位距离分布:")
        logger.info(f"    * 重叠 (0bp): {overlap_count}条")
        logger.info(f"    * 1-1000bp: {len(coloc_df[(coloc_df['distance'] > 0) & (coloc_df['distance'] <= 1000)])}条")
        logger.info(f"    * 1001-5000bp: {len(coloc_df[(coloc_df['distance'] > 1000) & (coloc_df['distance'] <= 5000)])}条")
        logger.info(f"    * 5001-10000bp: {len(coloc_df[(coloc_df['distance'] > 5000) & (coloc_df['distance'] <= 10000)])}条")
        
        # 最常见的组合
        if len(category_counts) > 0:
            top_category = category_counts.sort_values('count', ascending=False).head(5)
            logger.info("  - 最常见的ARG家族-MGE类别组合:")
            for _, row in top_category.iterrows():
                logger.info(f"    * {row['arg_family']} + {row['mge_category']}: {row['count']}次")
        
        # 最常见的亚类组合
        if len(subtype_counts) > 0:
            top_subtype = subtype_counts.sort_values('count', ascending=False).head(5)
            logger.info("  - 最常见的ARG亚类-MGE亚类组合:")
            for _, row in top_subtype.iterrows():
                logger.info(f"    * {row['arg_subtype']} + {row['mge_subtype']}: {row['count']}次")
        
        logger.info("共定位分析报告已保存到:")
        logger.info("  - ./analysis_results-HDS/co_localization_details.csv")
        logger.info("  - ./analysis_results-HDS/co_localization_summary_counts.csv")
        logger.info("  - ./analysis_results-HDS/co_localization_detailed_counts.csv")
        logger.info("  - ./analysis_results-HDS/co_localization_subtype_counts.csv")
    
    except Exception as e:
        logger.error(f"生成共定位报告时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())

# 主执行函数
def main():
    # 定义文件路径
    gff_file = "./analysis_results-HDS/predicted_genes.gff"
    arg_file = "./analysis_results-HDS/ARG_card_results.m8"
    mge_file = "./analysis_results-HDS/MGE_pfam_results.domtblout"
    arg_db = "./analysis_results-HDS/arg_results.db"
    mge_db = "./analysis_results-HDS/mge_results.db"
    association_db = "./analysis_results-HDS/arg_mge_associations.db"
    
    # 处理GFF文件
    gff_df, k141_to_coords, gff_id_to_coords = process_gff_file(gff_file)
    if gff_df.empty:
        logger.error("GFF处理失败，终止")
        return
    
    # 处理ARG数据
    logger.info("\n处理ARG数据...")
    arg_data = process_arg_data(arg_file, k141_to_coords, gff_id_to_coords)
    
    # 创建ARG数据库
    if not arg_data.empty:
        create_arg_database(arg_data, arg_db)
    
    # 处理MGE数据并创建数据库
    logger.info("\n处理MGE数据...")
    mge_success = create_mge_database(mge_file, k141_to_coords, gff_id_to_coords, mge_db)
    
    # 创建ARG-MGE关联分析
    success, pairs_df, common_contigs = False, pd.DataFrame(), set()
    if mge_success and os.path.exists(arg_db):
        logger.info("\n创建ARG-MGE关联分析...")
        success, pairs_df, common_contigs = create_arg_mge_association(arg_db, mge_db, association_db)
    
    # 如果有共定位contig，执行分类学分析
    tax_dict = {}
    if success and common_contigs:
        logger.info("\n执行分类学分析...")
        tax_dict = perform_taxonomy_analysis(common_contigs)
        
        # 如果获取了分类学信息，更新关联数据库
        if tax_dict:
            logger.info("更新分类学信息到关联数据库...")
            update_taxonomy_in_database(tax_dict, association_db)
    
    # 生成共定位分析报告
    logger.info("\n生成共定位分析报告...")
    generate_colocalization_reports(pairs_df, tax_dict)
    
    logger.info("\n数据处理、关联分析、分类学分析和共定位分析全部完成!")

if __name__ == "__main__":
    main()
EOF

# 添加执行权限
chmod +x ./scripts/process_annotations_and_colocalization.py

# 执行综合处理脚本
python ./scripts/process_annotations_and_colocalization.py