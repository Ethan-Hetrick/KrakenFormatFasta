import re
import argparse

def parse_conifer_output(conifer_list: list, compiled_pattern, rtl_threshold: float = None) -> int:
    try:
        # This will trigger if there is an unclassified read
        if conifer_list[0] == 'U':
            return 2
        # If read is classified and classified correctly
        elif conifer_list[0] == 'C' and compiled_pattern.search(conifer_list[1]).group(1) == conifer_list[2]:
            return 1
        # If read was classified but not classified correctly
        elif conifer_list[0] == 'C' and compiled_pattern.search(conifer_list[1]).group(1) != conifer_list[2]:
            return 0

    except (IndexError, AttributeError):
        print(f"Error processing line: {conifer_list}")

def process_conifer_output(confidence_score: float, rtl_threshold: float, conifer_output_path: str):
    classified_correctly = 0
    classified_incorrectly = 0
    total_reads = 0
    unclassified = 0

    # Regular expression pattern
    taxid_pattern = re.compile(r'taxid\|(\d+)\|')

    # Compile the pattern
    compiled_pattern = re.compile(taxid_pattern)

    with open(conifer_output_path, 'r') as conifer_output:
        for line in conifer_output:
            total_reads += 1
            conifer_list = line.split()
            result = parse_conifer_output(conifer_list, compiled_pattern)

            if (confidence_score is not None and float(conifer_list[-4]) < confidence_score) or (rtl_threshold is not None and float(conifer_list[-1]) < rtl_threshold):
                unclassified += 1
                continue
            elif result == 1:
                classified_correctly += 1
            elif result == 0:
                classified_incorrectly += 1
            elif result == 2:
                unclassified += 1

    return classified_correctly, classified_incorrectly, total_reads, unclassified

def print_results(file, confidence, rtl_threshold, classified_correctly, classified_incorrectly, total_reads, unclassified):

    precision = classified_correctly / (classified_correctly + classified_incorrectly) if total_reads != 0 else 0
    recall = classified_correctly / total_reads if total_reads != 0 else 0

    print("file,confidence_cutoff,rtl_threshold,total_reads,classified_correctly,classified_incorrectly,unclassified,precision,recall")
    print(f"{file},{confidence},{rtl_threshold},{total_reads},{classified_correctly},{classified_incorrectly},{unclassified},{precision},{recall}")

def main():
    parser = argparse.ArgumentParser(description="Process Conifer output.")
    parser.add_argument("-c", "--confidence", type=float, default=None, help="Confidence score threshold")
    parser.add_argument("-r", "--rtl_threshold", type=float, default=None, help="RTL score threshold")
    parser.add_argument("-i", "--input", required=True, help="Path to Conifer output file")
    args = parser.parse_args()

    classified_correctly, classified_incorrectly, total_reads, unclassified = process_conifer_output(args.confidence, args.rtl_threshold, args.input)
    print_results(args.input, args.confidence, args.rtl_threshold, classified_correctly, classified_incorrectly, total_reads, unclassified)

if __name__ == "__main__":
    main()