import sys, os
import subprocess
import re
import json
import argparse
from datetime import datetime

class ScreenOutput:
    contents = ''
    @staticmethod
    def add(c):
        ScreenOutput.contents += f"{c}\n"
        print(c)
    @staticmethod
    def writeScreenOutputToFile(filename):
        with open(filename, 'w') as f:
            f.write(ScreenOutput.contents)

def stripSeq(s):
    return ''.join(s.split('\n')[1:])

class PcrTemplate:
    def __init__(self, s):
        seq = s.lower()
        center = len(s) // 2
        self.upperCase = s.upper()
        self.lowerCase = s.lower()
        self.center5 = seq[center - 2:center + 3].upper()
        self.highlightCenter5 = seq[:center - 2] + self.center5 + seq[center + 3:]
        self.template = self.highlightCenter5
    def mute(self, pos, mutAllele):
        return PcrTemplate(self.template[:pos] + mutAllele + self.template[pos + 1:])

def complement(c):
    return {'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
            'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}.get(c, c)

def complementString(s):
    return ''.join([complement(c) for c in s])

def parseFastaFile(filename):
    with open(filename, 'r') as f:
        content = f.read()
    lines = content.strip().split('\n')
    header = lines[0].lstrip('>')
    sequence = ''.join(lines[1:])
    snp_match = re.search(r'\[([ACGT])/([ACGT])\]', sequence)
    if not snp_match:
        ScreenOutput.add("Error: No SNP found in the FASTA sequence.")
        sys.exit(1)
    wild_type = snp_match.group(1)
    mutation = snp_match.group(2)
    snp_pos = snp_match.start()
    clean_sequence = sequence.replace(snp_match.group(0), wild_type)
    output_dir = header.split()[0]
    return clean_sequence, wild_type, mutation, snp_pos, output_dir

def armsTable(tem, alt, adjacent):
    pair = ord(tem) + ord(alt)
    if pair in [ord('A') + ord('A'), ord('G') + ord('G')]:
        return 'A' if adjacent in ['A', 'C'] else 'G'
    elif pair in [ord('A') + ord('G'), ord('T') + ord('C'), ord('C') + ord('C')]:
        return 'C' if adjacent in ['A', 'C'] else 'T'
    elif pair == ord('A') + ord('C'):
        return {'A': 'G', 'G': 'A', 'C': 'C', 'T': 'T'}[adjacent]
    elif pair == ord('T') + ord('T'):
        return 'C' if adjacent in ['A', 'C'] else 'T'
    elif pair == ord('T') + ord('G'):
        return {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}[adjacent]
    return ''

def armsPcrTemplate(w, m, secondMutDistance):
    centerPos = len(m.template) // 2
    wtAllele = w.template[centerPos]
    wtAlleleCom = complement(wtAllele)
    mutAllele = m.template[centerPos]
    mutAlleleCom = complement(mutAllele)
    ajacentLeftAllele = w.template[centerPos - secondMutDistance]
    ajacentLeftAlleleCom = complement(ajacentLeftAllele)
    ajacentRightAllele = w.template[centerPos + secondMutDistance]
    ajacentRightAlleleCom = complement(ajacentRightAllele)
    wtTemForceLeftMissmatch = armsTable(wtAllele, mutAlleleCom, ajacentLeftAlleleCom)
    mutTemForceLeftMissmatch = armsTable(mutAllele, wtAlleleCom, ajacentLeftAlleleCom)
    wtTemForceRightMissmatch = armsTable(wtAlleleCom, mutAllele, ajacentRightAllele)
    mutTemForceRightMissmatch = armsTable(mutAlleleCom, wtAllele, ajacentRightAllele)
    return m.mute(centerPos - secondMutDistance, wtTemForceLeftMissmatch).mute(centerPos + secondMutDistance, mutTemForceRightMissmatch)

def primer3(parameters, parameterFileName, primerFileName):
    with open(parameterFileName, 'w') as f:
        f.write(parameters)
    os.system(f"primer3_core -output={primerFileName} -format_output -strict_tags {parameterFileName}")

def getMaskedTemplateFromFasta(sequence, snp_pos, output_dir, lengthPcrTemplate, nomask):
    halfLength = lengthPcrTemplate // 2
    start = max(0, snp_pos - halfLength)
    end = min(len(sequence), snp_pos + halfLength + 1)
    padded_seq = 'N' * max(0, halfLength - snp_pos) + sequence[start:end] + 'N' * max(0, (snp_pos + halfLength + 1) - len(sequence))
    centered_pos = len(padded_seq) // 2
    offset = centered_pos - (snp_pos - start)
    centered_seq = ('N' * offset + padded_seq[:-offset]) if offset > 0 else (padded_seq[-offset:] + 'N' * (-offset)) if offset < 0 else padded_seq
    pcrTemplate = PcrTemplate(centered_seq)
    templateFilename = f"{output_dir}/{lengthPcrTemplate}.txt"
    with open(templateFilename, 'w') as f:
        f.write(pcrTemplate.highlightCenter5)
    ScreenOutput.add(f"Genome sequence around SNP ({lengthPcrTemplate} bp) = {pcrTemplate.highlightCenter5}")
    ScreenOutput.add(f"Central 5 bp = {pcrTemplate.center5}")
    outputName = f"{output_dir}/masked.template.txt"
    with open(outputName, 'w') as f:
        f.write(pcrTemplate.template)
    ScreenOutput.add("*** nomask == True, therefore common SNPs are not masked as 'n' ***" if nomask else "*** nomask == False, but no masking done for FASTA input ***")
    return pcrTemplate

def designPrimers(filename):
    with open(filename) as f:
        sequence = f.read()
        halfLengthPlusOne = str(len(sequence) // 2 + 1)
    para_file = 'parameterOverride.json' if os.path.exists('parameterOverride.json') else 'parameterDefault.json'
    with open(para_file) as f:
        paraUser = json.load(f)
    paraUserStr = ''.join(f'{key}={value}\n' for key, value in sorted(paraUser.items()))
    for end_type in ['right', 'left']:
        parameters = f"""SEQUENCE_ID={filename}
SEQUENCE_TEMPLATE={sequence}
SEQUENCE_FORCE_{end_type.upper()}_END={halfLengthPlusOne}
{paraUserStr}="""
        parameterFileName = filename + f'.{end_type}.parameters.txt'
        primerFileName = f"{filename}.{end_type}.primer3.txt"
        primer3(parameters, parameterFileName, primerFileName)
        with open(primerFileName) as f:
            oligoFileContents = f.readlines()
        ScreenOutput.add(f"=========Force {end_type.capitalize()}==========")
        ScreenOutput.add(''.join(oligoFileContents[:7]))

def getAllArmsTemplates(snpID, wtAllele, mutAllele, pcrMaskedTemplate):
    position = len(pcrMaskedTemplate.template) // 2
    filename = f"{snpID}/masked.template"
    with open(filename + '.wild.txt', 'w') as f:
        f.write(pcrMaskedTemplate.template)
    ScreenOutput.add(f"WildTemplate = {pcrMaskedTemplate.center5}")
    pcrMaskedMutTemplate = pcrMaskedTemplate.mute(position, mutAllele)
    with open(filename + '.mutation.txt', 'w') as f:
        f.write(pcrMaskedMutTemplate.template)
    ScreenOutput.add(f"MutationTemplate = {pcrMaskedMutTemplate.center5}")
    for dist in [1, 2]:
        ScreenOutput.add(f"\nAdditional mutation at -{dist + 1} position")
        mut = armsPcrTemplate(pcrMaskedTemplate, pcrMaskedMutTemplate, dist)
        mut_filename = filename + f'.mutationTemplate.minus{dist + 1}.txt'
        with open(mut_filename, 'w') as f:
            f.write(mut.template)
        ScreenOutput.add(f"\nMutation template (-{dist + 1})")
        designPrimers(mut_filename)

class ABMPCRPrimerDesigner:
    def __init__(self, sequence):
        self.sequence = re.sub(r'[^ACGT]', '', sequence.upper())
    def validate_sequence(self):
        if not self.sequence:
            return False, "Error: Sequence is empty"
        if not re.match(r'^[ATCG]+$', self.sequence):
            return False, f"Error: Sequence contains invalid characters: {set(self.sequence) - set('ACGT')}"
        if len(self.sequence) < 18 or len(self.sequence) > 27:
            return False, "Error: Sequence length should be 18-27 bases"
        return True, ""
    def process_dna(self):
        valid, msg = self.validate_sequence()
        if not valid:
            return msg
        return self.process_last_base()
    def process_last_base(self):
        n = len(self.sequence)
        result = f"ARMS-PCR Primer Sequence: {self.sequence}\nTotal bases: {n}\n\n"
        result += f"{self.sequence[:-1]}\033[38;2;255;0;0m{self.sequence[-1]}\033[0m\n\n"
        result += self.apply_rule_p(n)
        return result
    def apply_rule_p(self, n):
        start = 6
        end = n - 7
        result = ""
        for i in range(end, start - 1, -1):
            modified_sequence = list(self.sequence)
            if modified_sequence[i] in 'CGTA':
                color = '\033[38;2;0;0;255m'
                replacement = {'C': 'T', 'G': 'A', 'T': 'A', 'A': 'T'}[modified_sequence[i]]
                modified_sequence[i] = f"{color}{replacement}\033[0m"
                result += f"Modification at position {i+1}: {''.join(modified_sequence)}\n"
        return result or "No modifications made\n"

class FileOutputHandler:
    def __init__(self, output_path=None):
        self.output_path = output_path
        self.file_ext = os.path.splitext(output_path)[1].lower() if output_path else None
    def save_to_file(self, content):
        if not self.output_path:
            return
        try:
            if self.file_ext == '.html':
                self._save_to_html(content)
            else:
                with open(self.output_path, 'w', encoding='utf-8') as f:
                    f.write(re.sub(r'\033\[[0-9;]+m', '', content))
                print(f"Results saved to text file: {self.output_path}")
        except Exception as e:
            print(f"Error saving file: {str(e)}")
    def _save_to_html(self, content):
        html_content = content.replace("\033[38;2;255;0;0m", "<span style='color:red;'>")
        html_content = html_content.replace("\033[38;2;0;0;255m", "<span style='color:blue;'>")
        html_content = html_content.replace("\033[0m", "</span>")
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        html = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>ABM-PCR Primer Design Results</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .red {{ color: red; }}
                .blue {{ color: blue; }}
                pre {{ white-space: pre-wrap; }}
                h2 {{ color: #333; border-bottom: 1px solid #ccc; padding-bottom: 5px; }}
                hr {{ border: 0; border-top: 1px solid #eee; margin: 20px 0; }}
            </style>
        </head>
        <body>
            <h1>ABM-PCR Primer Design Results</h1>
            <p>Generated on: {timestamp}</p>
            <p>A total of {len(primers)} primer sequences were found</p>
            {all_html_content}
        </body>
        </html>
        """
        with open(self.output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"Results saved to HTML file: {self.output_path}")

def parse_arguments():
    parser = argparse.ArgumentParser(description="ABM-PCR Primer Design Tool (Supports HTML Output)")
    parser.add_argument("--version", action="version", version="ABM-PCR Primer Designer 1.0")
    parser.add_argument("-f", "--fasta", type=str, required=True, help="Input FASTA file (required)")
    parser.add_argument("--length", type=int, default=1001, help="Length of PCR template (default: 1001)")
    parser.add_argument("--nomask", action="store_true", help="Do not mask common SNPs in template")
    parser.add_argument("--param", type=str, help="Custom parameters from JSON file (default: parameterDefault.json)")
    parser.add_argument("-o", "--output", type=str, help="Output HTML file path")
    return parser.parse_args()

def extract_primers(file_path):
    primer_neg2_right = primer_neg2_left = primer_neg3_right = primer_neg3_left = None
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        current_section = None
        for i, line in enumerate(lines):
            line = line.strip()
            if "Mutation template (-2)" in line:
                current_section = "-2"
            elif "Mutation template (-3)" in line:
                current_section = "-3"
            if current_section == "-2":
                if "Force Right" in line:
                    for j in range(i + 1, len(lines)):
                        if "RIGHT PRIMER" in lines[j]:
                            primer_neg2_right = lines[j].split()[-1]
                            break
                elif "Force Left" in line:
                    for j in range(i + 1, len(lines)):
                        if "LEFT PRIMER" in lines[j]:
                            primer_neg2_left = lines[j].split()[-1]
                            break
            if current_section == "-3":
                if "Force Right" in line:
                    for j in range(i + 1, len(lines)):
                        if "RIGHT PRIMER" in lines[j]:
                            primer_neg3_right = lines[j].split()[-1]
                            break
                elif "Force Left" in line:
                    for j in range(i + 1, len(lines)):
                        if "LEFT PRIMER" in lines[j]:
                            primer_neg3_left = lines[j].split()[-1]
                            break
    except FileNotFoundError:
        print(f"Error: File {file_path} not found")
    except Exception as e:
        print(f"Error processing file: {e}")
    return (primer_neg2_right, primer_neg2_left, primer_neg3_right, primer_neg3_left)

def highlight_primer(primer, index):
    if not primer:
        return primer
    if index == 0 or index == 1:  # Primer 1 和 Primer 2
        if len(primer) >= 2:
            return f"{primer[:-2]}<span class='blue'>{primer[-2]}</span><span class='red'>{primer[-1]}</span>"
        elif len(primer) >= 1:
            return f"{primer[:-1]}<span class='red'>{primer[-1]}</span>"
    elif index == 2 or index == 3:  # Primer 3 和 Primer 4
        if len(primer) >= 3:
            return f"{primer[:-3]}<span class='blue'>{primer[-3]}</span>{primer[-2]}<span class='red'>{primer[-1]}</span>"
        elif len(primer) >= 1:
            return f"{primer[:-1]}<span class='red'>{primer[-1]}</span>"
    return primer

if __name__ == '__main__':
    args = parse_arguments()
    if len(sys.argv) == 1 or '--help' in sys.argv:
        parser = parse_arguments()
        parser.print_help()
        sys.exit(0)
    lengthPcrTemplate = args.length
    nomask = args.nomask
    param_file = args.param
    fasta_file = args.fasta
    output_path = args.output
    if param_file:
        if os.path.exists(param_file):
            os.rename(param_file, 'parameterOverride.json')
        else:
            print(f"Error: Parameter file {param_file} not found")
            sys.exit(1)
    try:
        sequence, wtAllele, mutAllele, snp_pos, output_dir = parseFastaFile(fasta_file)
    except Exception as e:
        print(f"Error parsing FASTA file: {e}")
        sys.exit(1)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        ScreenOutput.add(f"Directory '{output_dir}' created to store output files")
    else:
        ScreenOutput.add(f"Directory '{output_dir}' exists and will store output files")
    snp_info = f"Wild type allele: {wtAllele}\nMutation allele: {mutAllele}\nSNP position: {snp_pos}"
    with open(f"{output_dir}/snp_info.txt", 'w') as f:
        f.write(snp_info)
    ScreenOutput.add('Details of your targeted SNP:\n' + snp_info)
    pcrMaskedTemplate = getMaskedTemplateFromFasta(sequence, snp_pos, output_dir, lengthPcrTemplate, nomask)
    getAllArmsTemplates(output_dir, wtAllele, mutAllele, pcrMaskedTemplate)
    ScreenOutput.writeScreenOutputToFile(f"{output_dir}/screenOutput.txt")
    primers = extract_primers(f"{output_dir}/screenOutput.txt")
    print(f"\nSuccessfully extracted {len(primers)} primer sequences")
    all_html_content = ""
    for i, primer in enumerate(primers):
        if primer:
            designer = ABMPCRPrimerDesigner(primer)
            result = designer.process_dna()
            print(result)
            html_content = result.replace("\033[38;2;255;0;0m", "<span class='red'>")
            html_content = html_content.replace("\033[38;2;0;0;255m", "<span class='blue'>")
            html_content = html_content.replace("\033[0m", "</span>")
            highlighted_primer = highlight_primer(primer, i)
            all_html_content += f"<h2>Primer {i + 1}: {highlighted_primer}</h2><pre>{html_content}</pre><hr>"
    if output_path and all_html_content:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        full_html = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>ABM-PCR Primer Design Results</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .red {{ color: red; font-weight: bold; }}
                .blue {{ color: blue; font-weight: bold; }}
                pre {{ white-space: pre-wrap; }}
                h2 {{ color: #333; border-bottom: 1px solid #ccc; padding-bottom: 5px; }}
                hr {{ border: 0; border-top: 1px solid #eee; margin: 20px 0; }}
                .primer {{ font-family: monospace; font-size: 1.1em; }}
            </style>
        </head>
        <body>
            <h1>ABM-PCR Primer Design Results</h1>
            <p>Generated on: {timestamp}</p>
            <p>A total of {len(primers)} primer sequences were found</p>
            {all_html_content}
        </body>
        </html>
        """
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(full_html)
            print(f"All results have been merged and saved to HTML file: {output_path}")
        except Exception as e:
            print(f"Error saving HTML file: {str(e)}")
            print("Please check file permissions or try a different output path")
    elif output_path and not all_html_content:
        print("Warning: No primer sequences found, HTML file not generated.")
    elif all_html_content and not output_path:
        print("Note: Primer analysis results generated but no output path specified, not saved to file.")
        print("Use --output parameter to specify HTML output file path.")
