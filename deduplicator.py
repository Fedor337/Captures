import sys
from collections import defaultdict

def process_fasta(input_file, output_file):

    sequences = defaultdict(list)
    current_header = None
    current_sequence = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                 
                    full_seq = ''.join(current_sequence)
                    sequences[full_seq].append(current_header)
                current_header = line
                current_sequence = []
            else:

                cleaned_line = ''.join(line.split())
                current_sequence.append(cleaned_line)
        
       
        if current_header is not None:
            full_seq = ''.join(current_sequence)
            sequences[full_seq].append(current_header)


    unique_count = 0
    with open(output_file, 'w') as out:
        for seq, headers in sequences.items():

            out.write(f"{headers[0]}\n")

            for i in range(0, len(seq), 80):
                out.write(f"{seq[i:i+80]}\n")
            unique_count += 1


    total_sequences = sum(len(headers) for headers in sequences.values())
    duplicates_removed = total_sequences - unique_count

    print(f"Processed {total_sequences} total sequences")
    print(f"Found {unique_count} unique sequences")
    print(f"Removed {duplicates_removed} duplicates")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.fa> <output.fa>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_fasta(input_file, output_file)
