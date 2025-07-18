#!/bin/bash
# 文件名: 12_visualization_and_statistics.sh
# 用途: 创建可视化和统计分析脚本并执行

cat > ./scripts/visualize_and_analyze_results.py << 'EOF'
#!/usr/bin/env python3
# 首先设置matplotlib后端为非交互式Agg后端，必须在导入pyplot之前设置
import matplotlib
matplotlib.use('Agg')  # 非交互式后端，解决Qt错误

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import sys
import logging
from scipy.stats import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests
import sqlite3
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from adjustText import adjust_text
from matplotlib_venn import venn2, venn3
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils.dataframe import dataframe_to_rows
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("./logs/visualization_analysis.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# 确保输出目录存在
os.makedirs("./logs", exist_ok=True)
os.makedirs("./analysis_results/figures", exist_ok=True)
os.makedirs("./analysis_results/excel_data", exist_ok=True)

logger.info("开始结果可视化与统计分析...")

# 设置所有图表的字体为 Times New Roman
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
plt.rcParams['mathtext.fontset'] = 'stix'  # 与Times New Roman兼容的数学字体

#################################################
# 数据加载函数
#################################################

def load_data():
    """加载分析所需的所有数据文件"""
    result_files = {
        'summary_counts': './analysis_results/co_localization_summary_counts.csv',
        'detailed_counts': './analysis_results/co_localization_detailed_counts.csv',
        'subtype_counts': './analysis_results/co_localization_subtype_counts.csv',
        'details': './analysis_results/co_localization_details.csv',
        'association_db': './analysis_results/arg_mge_associations.db',
        'taxonomy': './analysis_results/contig_taxonomy.csv',
        'taxonomy_with_coloc': './analysis_results/co_localization_with_taxonomy.csv'
    }
    
    # 检查核心文件是否存在
    core_files = ['summary_counts', 'detailed_counts', 'details']
    missing_core = []
    for name in core_files:
        if not os.path.exists(result_files[name]):
            missing_core.append(result_files[name])
            
    if missing_core:
        logger.error(f"以下核心文件不存在: {', '.join(missing_core)}")
        logger.error("请先运行共定位分析脚本生成这些文件")
        return None, None, None, None, None
    
    # 读取核心数据
    try:
        summary_df = pd.read_csv(result_files['summary_counts'])
        detailed_df = pd.read_csv(result_files['detailed_counts'])
        
        # 尝试读取详细数据
        try:
            details_df = pd.read_csv(result_files['details'])
        except:
            # 如果CSV读取失败，尝试从数据库读取
            conn = sqlite3.connect(result_files['association_db'])
            details_df = pd.read_sql(
                "SELECT * FROM arg_mge_associations WHERE distance <= 10000", 
                conn
            )
            conn.close()
        
        # 尝试读取亚类数据
        subtype_df = None
        if os.path.exists(result_files['subtype_counts']) and os.path.getsize(result_files['subtype_counts']) > 0:
            subtype_df = pd.read_csv(result_files['subtype_counts'])
            logger.info(f"- 亚类计数数据: {len(subtype_df)}行")
        else:
            logger.warning("亚类数据文件不存在或为空，部分可视化功能将不可用")
        
        # 尝试读取分类学数据
        taxonomy_df = None
        # 首先尝试读取with_taxonomy文件，因为它包含共定位事件
        if os.path.exists(result_files['taxonomy_with_coloc']) and os.path.getsize(result_files['taxonomy_with_coloc']) > 0:
            taxonomy_df = pd.read_csv(result_files['taxonomy_with_coloc'])
            logger.info(f"- 共定位事件分类学数据: {len(taxonomy_df)}行")
        # 如果没有，尝试读取常规分类学文件
        elif os.path.exists(result_files['taxonomy']) and os.path.getsize(result_files['taxonomy']) > 0:
            taxonomy_df = pd.read_csv(result_files['taxonomy'])
            logger.info(f"- contig分类学数据: {len(taxonomy_df)}行")
        else:
            logger.warning("分类学数据文件不存在或为空，分类学可视化功能将不可用")
        
        # 检查分类学信息是否在details_df中
        if taxonomy_df is None and 'phylum' in details_df.columns and 'genus' in details_df.columns:
            # 使用details_df中的分类学信息
            taxonomy_df = details_df[['contig', 'phylum', 'genus']].drop_duplicates()
            logger.info(f"- 从共定位详情提取的分类学数据: {len(taxonomy_df)}行")
        
        logger.info(f"成功加载数据:")
        logger.info(f"- 汇总计数数据: {len(summary_df)}行")
        logger.info(f"- 详细计数数据: {len(detailed_df)}行")
        logger.info(f"- 共定位详细信息: {len(details_df)}行")
        
        return summary_df, detailed_df, details_df, subtype_df, taxonomy_df
        
    except Exception as e:
        logger.error(f"加载数据时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return None, None, None, None, None

#################################################
# 数据导出到Excel函数
#################################################

def export_to_excel(data, filename, sheet_name='Data'):
    """将Pandas DataFrame导出到格式化的Excel文件"""
    try:
        excel_path = f"./analysis_results/excel_data/{filename}.xlsx"
        
        # 创建新的Excel文件
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = sheet_name
        
        # 将DataFrame数据写入Excel
        for r in dataframe_to_rows(data, index=False, header=True):
            ws.append(r)
        
        # 格式化标题行
        for cell in ws[1]:
            cell.font = Font(bold=True, name='Times New Roman')
            cell.alignment = Alignment(horizontal='center')
            cell.fill = PatternFill(start_color="E0E0E0", end_color="E0E0E0", fill_type="solid")
        
        # 设置列宽自适应
        for column in ws.columns:
            max_length = 0
            column_letter = column[0].column_letter
            for cell in column:
                try:
                    if len(str(cell.value)) > max_length:
                        max_length = len(str(cell.value))
                except:
                    pass
            adjusted_width = (max_length + 2) * 1.2
            ws.column_dimensions[column_letter].width = adjusted_width
        
        # 设置所有单元格字体
        for row in ws.iter_rows(min_row=2):
            for cell in row:
                cell.font = Font(name='Times New Roman')
        
        # 保存Excel文件
        wb.save(excel_path)
        logger.info(f"数据已导出到Excel文件: {excel_path}")
        return True
    except Exception as e:
        logger.error(f"导出Excel时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

#################################################
# 可视化函数
#################################################

def create_bubble_plot(summary_df):
    """创建ARG-MGE共定位气泡图"""
    logger.info("创建气泡图...")
    
    if summary_df.empty:
        logger.warning("数据为空，无法创建气泡图")
        return False
    
    try:
        # 获取ARG家族和MGE类别的数量
        arg_families = summary_df['arg_family'].unique()
        mge_categories = summary_df['mge_category'].unique()
        
        logger.info(f"气泡图将包含{len(arg_families)}个ARG家族和{len(mge_categories)}个MGE类别")
        
        # 创建数据透视表
        pivot_df = summary_df.pivot_table(
            index='arg_family', 
            columns='mge_category', 
            values='count',
            fill_value=0
        )
        
        # 导出数据到Excel
        export_to_excel(pivot_df.reset_index(), "arg_mge_bubble_data", "Colocalization_Counts")
        
        # 调整数据和图表大小
        if len(arg_families) > 20 or len(mge_categories) > 10:
            # 太多类别，只保留最常见的一些
            top_args = summary_df.groupby('arg_family')['count'].sum().nlargest(20).index
            top_mges = summary_df.groupby('mge_category')['count'].sum().nlargest(10).index
            
            # 过滤数据透视表
            pivot_df = pivot_df.loc[pivot_df.index.isin(top_args), pivot_df.columns.isin(top_mges)]
            logger.info(f"数据过滤后保留了top {len(pivot_df.index)} ARG家族和 {len(pivot_df.columns)} MGE类别")
        
        # 设置绘图样式
        fig_height = max(12, len(pivot_df.index) * 0.4)
        fig_width = max(14, len(pivot_df.columns) * 1.0)
        plt.figure(figsize=(fig_width, fig_height))
        sns.set(font_scale=1.2)
        sns.set_style("whitegrid")
        
        # 创建用于气泡大小的数组
        sizes = pivot_df.values.flatten()
        sizes_log = np.log1p(sizes)  # 对数变换
        sizes_scaled = sizes_log * 100  # 缩放到合适的大小
        
        # 创建网格坐标
        x, y = np.meshgrid(np.arange(pivot_df.shape[1]), np.arange(pivot_df.shape[0]))
        x = x.flatten()
        y = y.flatten()
        
        # 创建颜色映射 - 使用自定义颜色
        custom_cmap = LinearSegmentedColormap.from_list("custom_viridis", 
                                                       [(0, '#440154'), 
                                                        (0.5, '#21918c'), 
                                                        (1, '#fde725')])
        colors = custom_cmap(plt.Normalize()(sizes))
        
        # 绘制气泡图
        scatter = plt.scatter(x, y, s=sizes_scaled, c=colors, alpha=0.8, edgecolors='gray', linewidths=0.5)
        
        # 添加数值标签（仅添加到大于0的气泡）
        for i, (xi, yi, size, sz) in enumerate(zip(x, y, sizes, sizes_scaled)):
            if size > 0:
                plt.annotate(
                    f"{size:.0f}",
                    (xi, yi),
                    ha='center', va='center',
                    fontsize=8 if size < 10 else 9,
                    color='black' if sz < 250 else 'white'
                )
        
        # 设置轴标签和刻度
        plt.xticks(np.arange(len(pivot_df.columns)), pivot_df.columns, rotation=45, ha='right', fontname='Times New Roman')
        plt.yticks(np.arange(len(pivot_df.index)), pivot_df.index, fontname='Times New Roman')
        plt.xlabel('MGE Category', fontsize=14, fontname='Times New Roman')
        plt.ylabel('ARG Family', fontsize=14, fontname='Times New Roman')
        plt.title('ARG-MGE Co-localization Frequency', fontsize=16, fontname='Times New Roman')
        
        # 添加颜色条
        cbar = plt.colorbar(scatter)
        cbar.set_label('Co-localization Count', fontsize=12, fontname='Times New Roman')
        
        # 添加图例
        handles, labels = [], []
        for count in [1, 5, 10, 50]:
            handles.append(plt.scatter([], [], s=np.log1p(count)*100, color='gray', alpha=0.6))
            labels.append(f'{count}')
        plt.legend(handles, labels, title="Count", loc="upper right", bbox_to_anchor=(1.15, 0), 
                  prop={'family': 'Times New Roman'}, title_fontproperties={'family': 'Times New Roman'})
        
        # 调整布局并保存
        plt.tight_layout()
        plt.savefig('./analysis_results/figures/arg_mge_colocalization_bubble.png', dpi=300, bbox_inches='tight')
        plt.savefig('./analysis_results/figures/arg_mge_colocalization_bubble.pdf', bbox_inches='tight')
        plt.close()
        
        logger.info("气泡图已保存到 ./analysis_results/figures/arg_mge_colocalization_bubble.png/.pdf")
        return True
    
    except Exception as e:
        logger.error(f"创建气泡图时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_heatmap(summary_df):
    """创建ARG-MGE共定位热图"""
    logger.info("创建热图...")
    
    if summary_df.empty:
        logger.warning("数据为空，无法创建热图")
        return False
    
    try:
        # 创建数据透视表
        pivot_df = summary_df.pivot_table(
            index='arg_family', 
            columns='mge_category', 
            values='count',
            fill_value=0
        )
        
        # 导出数据到Excel
        export_to_excel(pivot_df.reset_index(), "arg_mge_heatmap_data", "Colocalization_Counts")
        
        # 调整数据和图表大小
        if len(pivot_df.index) > 20 or len(pivot_df.columns) > 10:
            # 太多类别，只保留最常见的一些
            top_args = summary_df.groupby('arg_family')['count'].sum().nlargest(20).index
            top_mges = summary_df.groupby('mge_category')['count'].sum().nlargest(10).index
            
            # 过滤数据透视表
            pivot_df = pivot_df.loc[pivot_df.index.isin(top_args), pivot_df.columns.isin(top_mges)]
        
        # 设置图表尺寸
        fig_height = max(12, len(pivot_df.index) * 0.4)
        fig_width = max(14, len(pivot_df.columns) * 1.0)
        
        plt.figure(figsize=(fig_width, fig_height))
        
        # 创建带标签的热图
        ax = sns.heatmap(
            pivot_df, 
            cmap='viridis', 
            annot=True,              # 显示数值
            fmt=".0f",               # 使用整数格式
            linewidths=0.5,          # 单元格边框
            cbar_kws={'label': 'Co-localization Count', 'shrink': 0.8}
        )
        
        # 设置字体为Times New Roman
        for text in ax.texts:
            text.set_fontname('Times New Roman')
        
        plt.title('ARG-MGE Co-localization Heatmap', fontsize=16, fontname='Times New Roman')
        plt.xlabel('MGE Category', fontsize=14, fontname='Times New Roman')
        plt.ylabel('ARG Family', fontsize=14, fontname='Times New Roman')
        
        # 设置刻度标签字体
        plt.xticks(rotation=45, ha='right', fontname='Times New Roman')
        plt.yticks(rotation=0, fontname='Times New Roman')
        
        # 调整颜色条标签字体
        cbar = ax.collections[0].colorbar
        cbar.ax.set_ylabel('Co-localization Count', fontname='Times New Roman')
        cbar.ax.tick_params(labelsize=10)
        for tick in cbar.ax.get_yticklabels():
            tick.set_fontname('Times New Roman')
        
        plt.tight_layout()
        plt.savefig('./analysis_results/figures/arg_mge_colocalization_heatmap.png', dpi=300, bbox_inches='tight')
        plt.savefig('./analysis_results/figures/arg_mge_colocalization_heatmap.pdf', bbox_inches='tight')
        plt.close()
        
        logger.info("热图已保存到 ./analysis_results/figures/arg_mge_colocalization_heatmap.png/.pdf")
        return True
    
    except Exception as e:
        logger.error(f"创建热图时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_subtype_heatmap(subtype_df):
    """创建ARG-MGE亚类共定位热图"""
    logger.info("创建亚类热图...")
    
    if subtype_df is None or subtype_df.empty:
        logger.warning("亚类数据为空，无法创建亚类热图")
        return False
    
    try:
        # 只选取前15个最常见的ARG和MGE亚类
        top_arg_subtypes = subtype_df.groupby('arg_subtype')['count'].sum().nlargest(15).index
        top_mge_subtypes = subtype_df.groupby('mge_subtype')['count'].sum().nlargest(15).index
        
        # 过滤数据
        filtered_df = subtype_df[
            subtype_df['arg_subtype'].isin(top_arg_subtypes) & 
            subtype_df['mge_subtype'].isin(top_mge_subtypes)
        ]
        
        # 创建数据透视表
        pivot_df = filtered_df.pivot_table(
            index='arg_subtype', 
            columns='mge_subtype', 
            values='count',
            fill_value=0
        )
        
        # 导出数据到Excel
        export_to_excel(pivot_df.reset_index(), "arg_mge_subtype_heatmap_data", "Subtype_Colocalization")
        
        # 设置图表尺寸
        fig_height = max(12, len(pivot_df.index) * 0.4)
        fig_width = max(15, len(pivot_df.columns) * 0.8)
        
        plt.figure(figsize=(fig_width, fig_height))
        
        # 创建带标签的热图
        ax = sns.heatmap(
            pivot_df, 
            cmap='viridis', 
            annot=True,
            fmt=".0f",
            linewidths=0.5,
            cbar_kws={'label': 'Co-localization Count', 'shrink': 0.8}
        )
        
        # 设置标签字体
        for text in ax.texts:
            text.set_fontname('Times New Roman')
            text.set_fontsize(8)  # 亚类名称可能较长，字体稍小
        
        plt.title('ARG-MGE Subtype Co-localization Heatmap', fontsize=16, fontname='Times New Roman')
        plt.xlabel('MGE Subtype', fontsize=14, fontname='Times New Roman')
        plt.ylabel('ARG Subtype', fontsize=14, fontname='Times New Roman')
        
        # 设置刻度标签字体和旋转
        plt.xticks(rotation=45, ha='right', fontname='Times New Roman', fontsize=9)
        plt.yticks(rotation=0, fontname='Times New Roman', fontsize=9)
        
        # 调整颜色条标签字体
        cbar = ax.collections[0].colorbar
        cbar.ax.set_ylabel('Co-localization Count', fontname='Times New Roman')
        for tick in cbar.ax.get_yticklabels():
            tick.set_fontname('Times New Roman')
        
        plt.tight_layout()
        plt.savefig('./analysis_results/figures/arg_mge_subtype_heatmap.png', dpi=300, bbox_inches='tight')
        plt.savefig('./analysis_results/figures/arg_mge_subtype_heatmap.pdf', bbox_inches='tight')
        plt.close()
        
        logger.info("亚类热图已保存到 ./analysis_results/figures/arg_mge_subtype_heatmap.png/.pdf")
        return True
    
    except Exception as e:
        logger.error(f"创建亚类热图时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_top_bar_charts(summary_df):
    """创建Top ARG-MGE共定位关系柱状图"""
    logger.info("创建Top10共定位关系柱状图...")
    
    if summary_df.empty:
        logger.warning("数据为空，无法创建柱状图")
        return False
    
    try:
        # 对结果进行排序
        sorted_df = summary_df.sort_values('count', ascending=False)
        
        # 准备导出到Excel的数据
        excel_data = {
            'Top_Combinations': sorted_df.head(10),
            'Top_ARG_Families': summary_df.groupby('arg_family')['count'].sum().nlargest(10).reset_index(),
            'Top_MGE_Categories': summary_df.groupby('mge_category')['count'].sum().nlargest(10).reset_index()
        }
        
        # 导出数据到Excel
        for sheet_name, data in excel_data.items():
            export_to_excel(data, f"top10_{sheet_name.lower()}", sheet_name)
        
        # 创建两个图表：Top10 ARG-MGE组合和Top10 ARG家族
        fig = plt.figure(figsize=(20, 10))
        gs = gridspec.GridSpec(1, 2, width_ratios=[1.5, 1])
        
        # 1. Top10 ARG-MGE组合
        ax1 = plt.subplot(gs[0])
        top_combinations = sorted_df.head(10)
        
        # 创建标签
        labels = [f"{row['arg_family']} - {row['mge_category']}" for _, row in top_combinations.iterrows()]
        
        # 创建渐变颜色
        colors = sns.color_palette("viridis", len(top_combinations))
        
        # 绘制水平条形图
        bars = ax1.barh(
            range(len(labels)),
            top_combinations['count'],
            color=colors,
            edgecolor='gray',
            linewidth=0.5,
            alpha=0.9
        )
        
        # 添加数值标签
        for i, bar in enumerate(bars):
            width = bar.get_width()
            ax1.text(
                width + (width * 0.02),  # 略微偏移
                bar.get_y() + bar.get_height()/2,
                f"{width:,}",
                va='center',
                fontname='Times New Roman'
            )
        
        # 设置标签和字体
        ax1.set_yticks(range(len(labels)))
        ax1.set_yticklabels(labels, fontname='Times New Roman')
        ax1.set_xlabel('Co-localization Count', fontsize=14, fontname='Times New Roman')
        ax1.set_title('Top10 ARG-MGE Co-localization Combinations', fontsize=16, fontname='Times New Roman')
        
        # 设置x轴刻度标签字体
        for tick in ax1.get_xticklabels():
            tick.set_fontname('Times New Roman')
        
        # 2. Top10 ARG家族总计数
        ax2 = plt.subplot(gs[1])
        top_arg_families = summary_df.groupby('arg_family')['count'].sum().nlargest(10)
        
        # 创建渐变颜色
        colors2 = sns.color_palette("magma", len(top_arg_families))
        
        # 绘制水平条形图
        bars2 = ax2.barh(
            range(len(top_arg_families.index)),
            top_arg_families.values,
            color=colors2,
            edgecolor='gray',
            linewidth=0.5,
            alpha=0.9
        )
        
        # 添加数值标签
        for i, bar in enumerate(bars2):
            width = bar.get_width()
            ax2.text(
                width + (width * 0.02),  # 略微偏移
                bar.get_y() + bar.get_height()/2,
                f"{width:,}",
                va='center',
                fontname='Times New Roman'
            )
        
        # 设置标签和字体
        ax2.set_yticks(range(len(top_arg_families.index)))
        ax2.set_yticklabels(top_arg_families.index, fontname='Times New Roman')
        ax2.set_xlabel('Total Co-localization Count', fontsize=14, fontname='Times New Roman')
        ax2.set_title('Top10 ARG Families (Total Count)', fontsize=16, fontname='Times New Roman')
        
        # 设置x轴刻度标签字体
        for tick in ax2.get_xticklabels():
            tick.set_fontname('Times New Roman')
            
        plt.tight_layout()
        plt.savefig('./analysis_results/figures/top10_colocalizations.png', dpi=300, bbox_inches='tight')
        plt.savefig('./analysis_results/figures/top10_colocalizations.pdf', bbox_inches='tight')
        plt.close()
        
        logger.info("Top10共定位关系柱状图已保存到 ./analysis_results/figures/top10_colocalizations.png/.pdf")
        
        # 创建第二个图：Top10 MGE类别总计数
        plt.figure(figsize=(12, 8))
        top_mge_categories = summary_df.groupby('mge_category')['count'].sum().nlargest(10)
        
        # 创建渐变颜色
        colors3 = sns.color_palette("crest", len(top_mge_categories))
        
        # 绘制条形图
        bars3 = plt.bar(
            top_mge_categories.index,
            top_mge_categories.values,
            color=colors3,
            edgecolor='gray',
            linewidth=0.5,
            alpha=0.9
        )
        
        # 添加数值标签
        for bar in bars3:
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width()/2,
                height + (height * 0.02),
                f"{height:,}",
                ha='center',
                fontname='Times New Roman'
            )
        
        # 设置字体
        plt.xticks(rotation=45, ha='right', fontname='Times New Roman')
        plt.yticks(fontname='Times New Roman')
        plt.ylabel('Total Co-localization Count', fontsize=14, fontname='Times New Roman')
        plt.title('Top10 MGE Categories (Total Count)', fontsize=16, fontname='Times New Roman')
        plt.tight_layout()
        plt.savefig('./analysis_results/figures/top10_mge_categories.png', dpi=300, bbox_inches='tight')
        plt.savefig('./analysis_results/figures/top10_mge_categories.pdf', bbox_inches='tight')
        plt.close()
        
        logger.info("Top10 MGE类别柱状图已保存到 ./analysis_results/figures/top10_mge_categories.png/.pdf")
        return True
    
    except Exception as e:
        logger.error(f"创建柱状图时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_distance_histogram(details_df):
    """创建距离分布直方图"""
    logger.info("创建共定位距离分布直方图...")
    
    if details_df.empty:
        logger.warning("数据为空，无法创建距离分布图")
        return False
    
    try:
        # 确保距离列存在
        if 'distance' not in details_df.columns:
            logger.error("详细数据中缺少'distance'列")
            return False
        
        # 定义距离区间
        bins = [0, 100, 500, 1000, 2000, 5000, 10000]
        bin_labels = ['0', '1-100', '101-500', '501-1000', '1001-2000', '2001-5000', '5001-10000']
        
        # 计算每个区间的计数
        counts = []
        for i in range(len(bins)-1):
            if i == 0:
                # 第一个区间只包含距离为0的（完全重叠）
                count = len(details_df[details_df['distance'] == 0])
            else:
                # 其他区间
                count = len(details_df[(details_df['distance'] > bins[i]) & (details_df['distance'] <= bins[i+1])])
            counts.append(count)
        
        # 确保bin_labels和counts长度匹配
        if len(bin_labels) != len(counts):
            bin_labels = bin_labels[:len(counts)]  # 截断bin_labels以匹配counts长度
        
        # 创建导出到Excel的数据
        distance_data = pd.DataFrame({
            'Distance_Range_bp': bin_labels,
            'Count': counts,
            'Percentage': [count/sum(counts)*100 for count in counts]
        })
        export_to_excel(distance_data, "distance_distribution", "Distance_Distribution")
        
        # 创建图形
        plt.figure(figsize=(14, 8))
        
        # 创建渐变颜色
        colors = sns.color_palette("viridis", len(counts))
        
        # 绘制条形图
        bars = plt.bar(bin_labels, counts, color=colors, edgecolor='black', linewidth=0.5)
        
        # 添加数值标签
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                plt.text(
                    bar.get_x() + bar.get_width()/2,
                    height + (max(counts) * 0.02),
                    f"{height:,}",
                    ha='center',
                    va='bottom',
                    fontname='Times New Roman'
                )
        
        # 设置轴标签和字体
        plt.xlabel('Distance Range (bp)', fontsize=14, fontname='Times New Roman')
        plt.ylabel('Co-localization Event Count', fontsize=14, fontname='Times New Roman')
        plt.title('ARG-MGE Co-localization Distance Distribution', fontsize=16, fontname='Times New Roman')
        
        # 设置刻度标签字体
        plt.xticks(fontname='Times New Roman')
        plt.yticks(fontname='Times New Roman')
        
        # 添加总计和百分比信息
        total = sum(counts)
        overlap_percent = counts[0] / total * 100 if total > 0 else 0
        
        plt.annotate(
            f"Total: {total:,} co-localization events\n"
            f"Complete overlap: {counts[0]:,} ({overlap_percent:.1f}%)",
            xy=(0.02, 0.95),
            xycoords='axes fraction',
            bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="gray", alpha=0.8),
            fontname='Times New Roman'
        )
        
        plt.tight_layout()
        plt.savefig('./analysis_results/figures/distance_distribution.png', dpi=300, bbox_inches='tight')
        plt.savefig('./analysis_results/figures/distance_distribution.pdf', bbox_inches='tight')
        plt.close()
        
        logger.info("距离分布直方图已保存到 ./analysis_results/figures/distance_distribution.png/.pdf")
        return True
    
    except Exception as e:
        logger.error(f"创建距离分布图时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_taxonomy_visualizations(details_df, taxonomy_df):
    """创建分类学相关的可视化图表"""
    logger.info("创建分类学相关可视化...")
    
    if taxonomy_df is None or taxonomy_df.empty:
        logger.warning("分类学数据为空，无法创建分类学相关可视化")
        return False
    
    try:
        # 1. 创建门级别分布饼图
        create_taxonomy_pie_chart(details_df, taxonomy_df, level='phylum')
        
        # 2. 创建属级别Top10条形图
        create_taxonomy_bar_chart(details_df, taxonomy_df, level='genus')
        
        # 3. 创建ARG家族-分类学热图
        create_arg_taxonomy_heatmap(details_df, taxonomy_df)
        
        return True
    
    except Exception as e:
        logger.error(f"创建分类学可视化时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_taxonomy_pie_chart(details_df, taxonomy_df, level='phylum'):
    """创建指定分类学级别的饼图"""
    logger.info(f"创建{level}级别分布饼图...")
    
    try:
        # 准备数据
        if 'phylum' in details_df.columns and 'genus' in details_df.columns:
            # 如果details_df已包含分类学信息，直接使用
            counts = details_df[level].value_counts()
        else:
            # 否则，合并分类学信息到details_df
            merged_df = pd.merge(
                details_df[['contig']], 
                taxonomy_df[['contig', level]],
                on='contig',
                how='left'
            )
            counts = merged_df[level].value_counts()
        
        # 处理缺失值
        if counts.index.isnull().any():
            counts = counts.fillna('unknown')
        
        # 处理'unknown'
        if 'unknown' in counts.index:
            unknown_count = counts['unknown']
            counts = counts.drop('unknown')
            counts['Unknown'] = unknown_count
        
        # 简化饼图：只显示Top10，其余归为"Others"
        if len(counts) > 10:
            top_taxa = counts.nlargest(9).index.tolist()
            others_sum = counts[~counts.index.isin(top_taxa)].sum()
            
            # 创建新的Series只包含Top9和Others
            new_counts = counts[counts.index.isin(top_taxa)]
            new_counts['Others'] = others_sum
            counts = new_counts
        
        # 导出数据到Excel
        export_data = counts.reset_index()
        export_data.columns = [level.capitalize(), 'Count']
        export_to_excel(export_data, f"{level}_distribution", f"{level.capitalize()}_Distribution")
        
        # 创建饼图
        plt.figure(figsize=(12, 8))
        
        # 使用鲜艳的颜色方案
        colors = plt.cm.tab20(np.linspace(0, 1, len(counts)))
        
        # 绘制饼图
        wedges, texts, autotexts = plt.pie(
            counts.values, 
            labels=counts.index, 
            autopct='%1.1f%%',
            textprops={'fontname': 'Times New Roman'},
            colors=colors,
            startangle=90,
            shadow=False,
            wedgeprops={'edgecolor': 'w', 'linewidth': 1},
            pctdistance=0.85
        )
        
        # 调整标签和百分比的字体
        for text in texts:
            text.set_fontname('Times New Roman')
        for autotext in autotexts:
            autotext.set_fontname('Times New Roman')
            autotext.set_fontsize(9)
        
        # 添加圆环效果
        centre_circle = plt.Circle((0, 0), 0.5, fc='white', edgecolor='gray')
        plt.gca().add_artist(centre_circle)
        
        # 添加标题
        plt.title(f'Distribution of Co-localization Events by {level.capitalize()}', 
                 fontsize=16, fontname='Times New Roman')
        
        # 添加注释展示总数
        total = counts.sum()
        plt.annotate(
            f"Total: {total:,} events",
            xy=(0.5, 0.5),
            xycoords='figure fraction',
            ha='center',
            fontsize=12,
            fontname='Times New Roman'
        )
        
        plt.tight_layout()
        plt.savefig(f'./analysis_results/figures/{level}_distribution_pie.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'./analysis_results/figures/{level}_distribution_pie.pdf', bbox_inches='tight')
        plt.close()
        
        logger.info(f"{level}级别分布饼图已保存到 ./analysis_results/figures/{level}_distribution_pie.png/.pdf")
        return True
    
    except Exception as e:
        logger.error(f"创建{level}级别分布饼图时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_taxonomy_bar_chart(details_df, taxonomy_df, level='genus'):
    """创建指定分类学级别的Top10条形图"""
    logger.info(f"创建Top10 {level}级别条形图...")
    
    try:
        # 准备数据
        if 'phylum' in details_df.columns and 'genus' in details_df.columns:
            # 如果details_df已包含分类学信息，直接使用
            counts = details_df[level].value_counts()
        else:
            # 否则，合并分类学信息到details_df
            merged_df = pd.merge(
                details_df[['contig']], 
                taxonomy_df[['contig', level]],
                on='contig',
                how='left'
            )
            counts = merged_df[level].value_counts()
        
        # 处理缺失值和unknown
        if counts.index.isnull().any():
            counts = counts.fillna('unknown')
        
        if 'unknown' in counts.index:
            counts = counts.rename({'unknown': 'Unknown'})
        
        # 获取Top10（不包括Unknown）
        if 'Unknown' in counts.index:
            unknown_count = counts['Unknown']
            counts = counts.drop('Unknown')
            top_taxa = counts.nlargest(9).index.tolist()
            top_counts = counts[counts.index.isin(top_taxa)]
            top_counts['Unknown'] = unknown_count
        else:
            top_taxa = counts.nlargest(10).index.tolist()
            top_counts = counts[counts.index.isin(top_taxa)]
        
        # 导出数据到Excel
        export_data = top_counts.reset_index()
        export_data.columns = [level.capitalize(), 'Count']
        export_to_excel(export_data, f"top10_{level}_distribution", f"Top10_{level.capitalize()}")
        
        # 创建条形图
        plt.figure(figsize=(14, 8))
        
        # 创建渐变颜色
        colors = plt.cm.tab20(np.linspace(0, 1, len(top_counts)))
        
        # 绘制条形图
        bars = plt.bar(
            top_counts.index,
            top_counts.values,
            color=colors,
            edgecolor='black',
            linewidth=0.5
        )
        
        # 添加数值标签
        for bar in bars:
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width()/2,
                height + (max(top_counts.values) * 0.02),
                f"{height:,}",
                ha='center',
                va='bottom',
                fontname='Times New Roman'
            )
        
        # 设置轴标签和字体
        plt.xlabel(level.capitalize(), fontsize=14, fontname='Times New Roman')
        plt.ylabel('Co-localization Event Count', fontsize=14, fontname='Times New Roman')
        plt.title(f'Top {len(top_counts)} {level.capitalize()} in ARG-MGE Co-localization Events', 
                 fontsize=16, fontname='Times New Roman')
        
        # 设置刻度标签字体并旋转
        plt.xticks(rotation=45, ha='right', fontname='Times New Roman')
        plt.yticks(fontname='Times New Roman')
        
        plt.tight_layout()
        plt.savefig(f'./analysis_results/figures/top10_{level}_distribution_bar.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'./analysis_results/figures/top10_{level}_distribution_bar.pdf', bbox_inches='tight')
        plt.close()
        
        logger.info(f"Top10 {level}级别条形图已保存到 ./analysis_results/figures/top10_{level}_distribution_bar.png/.pdf")
        return True
    
    except Exception as e:
        logger.error(f"创建Top10 {level}级别条形图时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_arg_taxonomy_heatmap(details_df, taxonomy_df):
    """创建ARG家族与分类学分布的热图"""
    logger.info("创建ARG家族-分类学热图...")
    
    try:
        # 准备数据
        if 'arg_family' in details_df.columns and 'phylum' in details_df.columns:
            # 如果details_df已包含分类学信息，直接使用
            df = details_df
        else:
            # 否则，合并分类学信息到details_df
            df = pd.merge(
                details_df, 
                taxonomy_df[['contig', 'phylum']],
                on='contig',
                how='left'
            )
        
        # 过滤掉未知门
        df = df.copy()
        df.loc[df['phylum'].isnull(), 'phylum'] = 'Unknown'
        df.loc[df['phylum'] == 'unknown', 'phylum'] = 'Unknown'
        
        # 获取Top ARG家族
        top_arg_families = df['arg_family'].value_counts().nlargest(15).index.tolist()
        
        # 获取Top Phyla（门）
        top_phyla = df['phylum'].value_counts().nlargest(10).index.tolist()
        if 'Unknown' in top_phyla:
            top_phyla.remove('Unknown')
            top_phyla.append('Unknown')  # 将Unknown放在最后
        
        # 过滤数据
        filtered_df = df[df['arg_family'].isin(top_arg_families) & df['phylum'].isin(top_phyla)]
        
        # 创建计数矩阵
        pivot_df = filtered_df.groupby(['arg_family', 'phylum']).size().unstack(fill_value=0)
        
        # 按总数重新排序ARG家族
        arg_family_order = filtered_df['arg_family'].value_counts().index.tolist()
        pivot_df = pivot_df.reindex(arg_family_order)
        
        # 导出数据到Excel
        export_to_excel(pivot_df.reset_index(), "arg_phylum_heatmap_data", "ARG_Phylum_Distribution")
        
        # 创建热图
        plt.figure(figsize=(15, 10))
        
        # 绘制热图
        ax = sns.heatmap(
            pivot_df, 
            cmap='viridis', 
            annot=True,
            fmt="d",
            linewidths=0.5,
            cbar_kws={'label': 'Co-localization Count', 'shrink': 0.8}
        )
        
        # 设置字体为Times New Roman
        for text in ax.texts:
            text.set_fontname('Times New Roman')
        
        plt.title('ARG Family Distribution Across Phyla', fontsize=16, fontname='Times New Roman')
        plt.xlabel('Phylum', fontsize=14, fontname='Times New Roman')
        plt.ylabel('ARG Family', fontsize=14, fontname='Times New Roman')
        
        # 设置刻度标签字体
        plt.xticks(rotation=45, ha='right', fontname='Times New Roman')
        plt.yticks(fontname='Times New Roman')
        
        # 调整颜色条标签字体
        cbar = ax.collections[0].colorbar
        cbar.ax.set_ylabel('Co-localization Count', fontname='Times New Roman')
        for tick in cbar.ax.get_yticklabels():
            tick.set_fontname('Times New Roman')
        
        plt.tight_layout()
        plt.savefig('./analysis_results/figures/arg_phylum_heatmap.png', dpi=300, bbox_inches='tight')
        plt.savefig('./analysis_results/figures/arg_phylum_heatmap.pdf', bbox_inches='tight')
        plt.close()
        
        logger.info("ARG家族-分类学热图已保存到 ./analysis_results/figures/arg_phylum_heatmap.png/.pdf")
        return True
    
    except Exception as e:
        logger.error(f"创建ARG家族-分类学热图时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

#################################################
# 统计分析函数
#################################################

def perform_fisher_test(summary_df, output_file='./analysis_results/fisher_test_results.csv', 
                      significant_file='./analysis_results/significant_associations.csv'):
    """执行Fisher精确检验分析ARG-MGE关联的显著性"""
    logger.info("开始执行Fisher精确检验...")
    
    if summary_df.empty:
        logger.warning("数据为空，无法进行Fisher检验")
        return False
    
    try:
        # 计算总共定位事件数
        total_coloc = summary_df['count'].sum()
        logger.info(f"共定位事件总数: {total_coloc}")
        
        if total_coloc == 0:
            logger.warning("没有检测到共定位事件，无法进行统计分析")
            return False
        
        # 准备Fisher精确检验结果
        fisher_results = []
        
        # 获取唯一的ARG家族和MGE类别
        arg_families = summary_df['arg_family'].unique()
        mge_categories = summary_df['mge_category'].unique()
        
        logger.info(f"分析{len(arg_families)}个ARG家族和{len(mge_categories)}个MGE类别间的关系...")
        
        # 对每个ARG-MGE组合进行测试
        for arg_family in arg_families:
            for mge_category in mge_categories:
                # 构建2x2列联表
                a = summary_df[(summary_df['arg_family'] == arg_family) & 
                            (summary_df['mge_category'] == mge_category)]['count'].sum()
                b = summary_df[summary_df['arg_family'] == arg_family]['count'].sum() - a
                c = summary_df[summary_df['mge_category'] == mge_category]['count'].sum() - a
                d = total_coloc - a - b - c
                
                # 检查是否有共定位事件
                if a > 0:
                    # 计算Fisher精确检验
                    table = np.array([[a, b], [c, d]])
                    odds_ratio, p_value = fisher_exact(table, alternative='greater')
                    
                    # 计算期望值 (行总和 * 列总和 / 总数)
                    expected = (a+b)*(a+c)/total_coloc if total_coloc > 0 else 0
                    
                    # 保存结果
                    fisher_results.append({
                        'arg_family': arg_family,
                        'mge_category': mge_category,
                        'observed': a,
                        'expected': expected,
                        'enrichment': a / expected if expected > 0 else float('inf'),
                        'odds_ratio': odds_ratio,
                        'p_value': p_value
                    })
        
        # 转换为DataFrame
        fisher_df = pd.DataFrame(fisher_results)
        
        # 多重检验校正
        if len(fisher_df) > 0:
            _, fisher_df['p_adjusted'], _, _ = multipletests(fisher_df['p_value'], method='fdr_bh')
            
            # 显示前几行
            logger.info("Fisher检验结果示例（前5行）:")
            for i, row in fisher_df.head(5).iterrows():
                logger.info(f"  {row['arg_family']} - {row['mge_category']}: "
                          f"观察值={row['observed']:.0f}, 期望值={row['expected']:.2f}, "
                          f"富集倍数={row['enrichment']:.2f}, "
                          f"调整后p值={row['p_adjusted']:.4f}")
                          
            # 保存完整结果
            fisher_df.to_csv(output_file, index=False)
            logger.info(f"全部Fisher检验结果已保存到: {output_file}")
            
            # 保存Excel格式
            export_to_excel(fisher_df, "fisher_test_results", "Fisher_Test_Results")
            
            # 筛选显著结果
            significant = fisher_df[fisher_df['p_adjusted'] < 0.05]
            significant = significant.sort_values('odds_ratio', ascending=False)
            
            if len(significant) > 0:
                significant.to_csv(significant_file, index=False)
                logger.info(f"显著关联结果已保存到: {significant_file}")
                
                # 保存Excel格式
                export_to_excel(significant, "significant_associations", "Significant_Associations")
                
                # 显示前几个显著结果
                logger.info(f"发现{len(significant)}个统计显著的关联，前5个显著关联:")
                for i, row in significant.head(5).iterrows():
                    logger.info(f"  {row['arg_family']} - {row['mge_category']}: "
                              f"观察值={row['observed']:.0f}, 期望值={row['expected']:.2f}, "
                              f"富集倍数={row['enrichment']:.2f}, "
                              f"调整后p值={row['p_adjusted']:.6f}")
                
                # 创建显著关联的可视化
                create_significant_association_plot(significant)
                
                return True
            else:
                logger.info("未发现统计显著的关联")
                return True
        else:
            logger.warning("没有足够的数据进行Fisher精确检验")
            return False
    
    except Exception as e:
        logger.error(f"执行Fisher检验时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_significant_association_plot(significant_df):
    """为显著关联创建可视化"""
    logger.info("创建显著关联图...")
    
    if significant_df.empty:
        logger.warning("没有显著关联，无法创建图表")
        return False
    
    try:
        # 只取前20个最显著的结果（如果超过20个）
        if len(significant_df) > 20:
            plot_df = significant_df.head(20)
        else:
            plot_df = significant_df
        
        # 创建组合标签
        plot_df['label'] = plot_df.apply(lambda x: f"{x['arg_family']} - {x['mge_category']}", axis=1)
        
        # 按照富集倍数排序
        plot_df = plot_df.sort_values('enrichment')
        
        # 创建水平条形图
        fig, ax = plt.subplots(figsize=(14, max(8, len(plot_df) * 0.4)))
        
        # 创建颜色映射基于p值
        colors = plt.cm.viridis(plt.Normalize()(plot_df['p_adjusted']))
        
        # 绘制条形图
        bars = ax.barh(plot_df['label'], plot_df['enrichment'], color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
        
        # 添加数值标签
        for i, bar in enumerate(bars):
            width = bar.get_width()
            ax.text(
                width + 0.1,
                bar.get_y() + bar.get_height()/2,
                f"{width:.2f}x",
                va='center',
                fontname='Times New Roman'
            )
        
        # 添加颜色条 - 基于p值的显著性
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=0, vmax=0.05))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('Adjusted p-value (FDR)', fontname='Times New Roman')
        
        # 设置刻度标签字体
        for tick in cbar.ax.get_yticklabels():
            tick.set_fontname('Times New Roman')
        
        # 设置轴标签和字体
        ax.set_xlabel('Enrichment Factor', fontsize=14, fontname='Times New Roman')
        ax.set_title('Significant ARG-MGE Associations (p < 0.05)', fontsize=16, fontname='Times New Roman')
        ax.grid(axis='x', linestyle='--', alpha=0.7)
        
        # 设置y轴刻度标签字体
        for tick in ax.get_yticklabels():
            tick.set_fontname('Times New Roman')
        
        # 设置x轴刻度标签字体
        for tick in ax.get_xticklabels():
            tick.set_fontname('Times New Roman')
        
        plt.tight_layout()
        plt.savefig('./analysis_results/figures/significant_associations.png', dpi=300, bbox_inches='tight')
        plt.savefig('./analysis_results/figures/significant_associations.pdf', bbox_inches='tight')
        plt.close()
        
        logger.info("显著关联图已保存到 ./analysis_results/figures/significant_associations.png/.pdf")
        return True
    
    except Exception as e:
        logger.error(f"创建显著关联图时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_network_visualization(details_df, subtype_df):
    """创建ARG-MGE网络可视化图"""
    logger.info("创建ARG-MGE网络图...")
    
    if details_df.empty:
        logger.warning("数据为空，无法创建网络图")
        return False
    
    # 网络图需要使用networkx库，但不是必须的，所以在这里检查是否有安装
    try:
        import networkx as nx
    except ImportError:
        logger.warning("未安装networkx库，无法创建网络图。可以使用pip install networkx安装")
        return False
    
    try:
        # 准备ARG-MGE链接数据
        if subtype_df is not None and not subtype_df.empty:
            # 使用亚类数据，更详细的网络图
            network_df = subtype_df[['arg_subtype', 'mge_subtype', 'count']].rename(
                columns={'arg_subtype': 'source', 'mge_subtype': 'target', 'count': 'weight'}
            )
            # 只保留权重最高的200条边，避免图太复杂
            if len(network_df) > 200:
                network_df = network_df.nlargest(200, 'weight')
            title = 'ARG-MGE Subtype Network'
        else:
            # 使用家族数据
            network_df = details_df.groupby(['arg_family', 'mge_category']).size().reset_index()
            network_df.columns = ['source', 'target', 'weight']
            # 只保留权重最高的100条边，避免图太复杂
            if len(network_df) > 100:
                network_df = network_df.nlargest(100, 'weight')
            title = 'ARG-MGE Family Network'
        
        # 导出网络数据到Excel
        export_to_excel(network_df, "arg_mge_network_data", "Network_Edges")
        
        # 创建有向图
        G = nx.DiGraph()
        
        # 添加节点
        arg_nodes = network_df['source'].unique()
        mge_nodes = network_df['target'].unique()
        
        # 添加ARG节点
        for node in arg_nodes:
            G.add_node(node, type='ARG')
        
        # 添加MGE节点
        for node in mge_nodes:
            G.add_node(node, type='MGE')
        
        # 添加边
        for _, row in network_df.iterrows():
            G.add_edge(row['source'], row['target'], weight=row['weight'])
        
        # 创建图
        plt.figure(figsize=(18, 12))
        
        # 使用spring布局
        pos = nx.spring_layout(G, k=0.15, iterations=50, seed=42)
        
        # 节点大小基于度
        node_degrees = dict(G.degree())
        node_sizes = [v * 30 for v in node_degrees.values()]
        
        # 根据节点类型设置颜色
        node_colors = []
        for node in G.nodes():
            if G.nodes[node]['type'] == 'ARG':
                node_colors.append('#ff7f0e')  # 橙色 for ARG
            else:
                node_colors.append('#1f77b4')  # 蓝色 for MGE
        
        # 绘制节点
        nx.draw_networkx_nodes(
            G, pos, 
            node_size=node_sizes,
            node_color=node_colors,
            alpha=0.8
        )
        
        # 根据权重设置边的宽度和透明度
        edges = G.edges()
        edge_weights = [G[u][v]['weight'] for u, v in edges]
        max_weight = max(edge_weights)
        edge_widths = [1.0 + 3.0 * (weight / max_weight) for weight in edge_weights]
        edge_alphas = [0.3 + 0.7 * (weight / max_weight) for weight in edge_weights]
        
        # 绘制边
        nx.draw_networkx_edges(
            G, pos, 
            width=edge_widths,
            alpha=0.5,
            edge_color='gray',
            arrows=False
        )
        
        # 为节点添加标签，只标注度大的节点
        node_labels = {}
        for node in G.nodes():
            if node_degrees[node] > 2:
                node_labels[node] = node
        
        # 绘制标签
        nx.draw_networkx_labels(
            G, pos, 
            labels=node_labels,
            font_size=10,
            font_family='Times New Roman',
            font_weight='bold'
        )
        
        # 添加图例
        arg_patch = Patch(color='#ff7f0e', label='ARG')
        mge_patch = Patch(color='#1f77b4', label='MGE')
        plt.legend(handles=[arg_patch, mge_patch], prop={'family': 'Times New Roman', 'size': 12})
        
        # 添加标题
        plt.title(title, fontsize=16, fontname='Times New Roman')
        
        # 关闭坐标轴
        plt.axis('off')
        
        plt.tight_layout()
        plt.savefig('./analysis_results/figures/arg_mge_network.png', dpi=300, bbox_inches='tight')
        plt.savefig('./analysis_results/figures/arg_mge_network.pdf', bbox_inches='tight')
        plt.close()
        
        logger.info("ARG-MGE网络图已保存到 ./analysis_results/figures/arg_mge_network.png/.pdf")
        return True
    
    except Exception as e:
        logger.error(f"创建网络图时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

#################################################
# 主函数
#################################################

def main():
    """主执行函数"""
    try:
        # 读取数据
        summary_df, detailed_df, details_df, subtype_df, taxonomy_df = load_data()
        if summary_df is None or detailed_df is None or details_df is None:
            logger.error("核心数据加载失败，无法继续分析")
            return
        
        # 执行可视化
        logger.info("\n===== 开始可视化结果 =====")
        
        # 创建气泡图
        create_bubble_plot(summary_df)
        
        # 创建热图
        create_heatmap(summary_df)
        
        # 创建亚类热图（如果有亚类数据）
        if subtype_df is not None and not subtype_df.empty:
            create_subtype_heatmap(subtype_df)
        
        # 创建柱状图
        create_top_bar_charts(summary_df)
        
        # 创建距离分布直方图
        create_distance_histogram(details_df)
        
        # 创建分类学可视化（如果有分类学数据）
        if taxonomy_df is not None and not taxonomy_df.empty:
            create_taxonomy_visualizations(details_df, taxonomy_df)
        
        # 创建网络图
        create_network_visualization(details_df, subtype_df)
        
        # 执行统计分析
        logger.info("\n===== 开始统计分析 =====")
        perform_fisher_test(summary_df)
        
        logger.info("\n===== 所有分析和可视化任务完成 =====")
    
    except Exception as e:
        logger.error(f"执行可视化和分析时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())

if __name__ == "__main__":
    main()
EOF

# 添加执行权限
chmod +x ./scripts/visualize_and_analyze_results.py

# 执行可视化和统计分析脚本
python ./scripts/visualize_and_analyze_results.py

echo "可视化与统计分析完成!"
echo "结果图像已保存到 ./analysis_results/figures/ 目录"
echo "图表数据已保存到 ./analysis_results/excel_data/ 目录"
echo "统计分析结果已保存到 ./analysis_results/ 目录"