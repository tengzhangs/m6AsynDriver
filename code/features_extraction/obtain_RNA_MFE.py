import os
import RNA
import pandas as pd
from pathlib import Path
from Bio import SeqIO

def validate_dot_bracket(sequence, structure):
    """
    验证点括号结构是否正确
    - 检查长度是否匹配
    - 检查括号是否配对
    """
    # 检查长度
    if len(sequence) != len(structure):
        return False, f"长度不匹配: 序列长度 {len(sequence)}, 结构长度 {len(structure)}"
    
    # 检查括号配对
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if not stack:
                return False, "括号不匹配: 多余的右括号"
            stack.pop()
    
    if stack:
        return False, "括号不匹配: 多余的左括号"
    
    return True, "结构有效"

def dna_to_rna(sequence):
    """将DNA序列转换为RNA序列（T -> U）"""
    if not isinstance(sequence, str):
        sequence = str(sequence)
    return sequence.upper().replace('T', 'U')

def process_rna_files(input_folder, output_folder):
    """处理指定文件夹中的所有FASTA文件"""
    
    print(f"开始处理文件夹: {input_folder}")
    
    # 检查输入文件夹是否存在
    if not os.path.exists(input_folder):
        print(f"错误: 输入文件夹 '{input_folder}' 不存在")
        return
    
    # 创建输出文件夹（如果不存在）
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    # 获取输入文件夹中的所有FASTA文件
    fasta_files = [file for file in os.listdir(input_folder) 
                  if file.endswith('.fa') or file.endswith('.fasta')]
    
    print(f"找到 {len(fasta_files)} 个FASTA文件")
    
    if len(fasta_files) == 0:
        print("没有找到FASTA文件，处理结束")
        return
    
    # 处理每个文件
    for file_name in fasta_files:
        input_path = os.path.join(input_folder, file_name)
        output_path = os.path.join(output_folder, os.path.splitext(file_name)[0] + '.csv')
        
        print(f"\n正在处理文件: {file_name}")
        
        try:
            # 读取FASTA文件
            records = list(SeqIO.parse(input_path, "fasta"))
            print(f"读取了 {len(records)} 条序列记录")
            
            # 初始化结果列表
            sequence_list = []
            structures = []
            energys = []
            # 处理每个序列
            for i, record in enumerate(records):
                seq_id = record.id
                original_seq = str(record.seq)
                sequence_list.append(original_seq)  # 保存原序列
                
                # 将DNA序列转换为RNA序列（T -> U）
                rna_sequence = dna_to_rna(original_seq)
                
                try:
                    # 使用ViennaRNA计算结构
                    structure, energy = RNA.fold(rna_sequence)
                    
                    # 验证结构
                    is_valid, message = validate_dot_bracket(rna_sequence, structure)
                    
                    # 打印结果
                    print(f"\n序列 {i+1}: {seq_id}")
                    print(f"序列: {original_seq}")
                    print(f"二级结构: {structure}")
                    print(f"自由能: {energy:.2f} kcal/mol")
                    print(f"结构有效: {'是' if is_valid else '否 - ' + message}")
                    
                    structures.append(structure)
                    
                    energys.append(energy)
                    
                    # 每处理10个序列显示一次进度
                    if (i + 1) % 10 == 0 and i + 1 < len(records):
                        print(f"已处理 {i+1}/{len(records)} 条序列")
                        
                except Exception as e:
                    print(f"处理序列 {seq_id} 时出错: {str(e)}")
                    # 添加默认结构
                    structures.append("." * len(rna_sequence))
            
            # 创建新的DataFrame
            result_df = pd.DataFrame({
                'x': sequence_list,  # 原始序列
                'y': structures,       # 点括号结构
                'MFE': energys   # 最小自由能
            })
            
            # 保存结果到新文件
            result_df.to_csv(output_path, index=False)
            print(f"已处理完成: {file_name} -> {output_path}")
            
        except Exception as e:
            print(f"处理文件 {file_name} 时出错: {str(e)}")

if __name__ == "__main__":
    try:
        # 测试RNA模块
        test_fold = RNA.fold("AUCG")
        print("ViennaRNA库初始化成功")
    except Exception as e:
        print(f"ViennaRNA初始化失败: {str(e)}")
        exit(1)
    
    # 保持原始的输入和输出文件夹路径
    #input_folder = r'D:\\research\\m6Acancer_prediction\\m6A_mutation\\mutation_data\\new_data\\driver_passenger\\RNA_fold\\pos'
    #input_folder = r'D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\RNA_fold\\Mu_VS_WT\\fast_file\\muta\\'
    input_folder = r'D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\RNA_fold\\Mu_VS_WT\\fast_file\\WT\\'
    output_folder = r'D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\RNA_fold\\Mu_VS_WT\\processed'
    #
    #input_folder = r'D:\\research\\m6Atranslation\data\\single_base_level\\cdhit_proces\\neg_RNA_fold'
    #output_folder = r'D:\\research\\m6Atranslation\data\\single_base_level\\cdhit_proces\\neg_RNA_fold\\processed'
    
    # 处理文件
    process_rna_files(input_folder, output_folder)
    print("\n处理完成!")
