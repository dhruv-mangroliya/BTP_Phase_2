import random
from collections import defaultdict

# Function to generate a random protein sequence of given length
def generate_protein_sequence(length):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # 20 standard amino acids
    return ''.join(random.choices(amino_acids, k=length))

# Function to fragment the protein sequence into chunks of max length 1000
def fragment_protein_sequence(sequence, max_length=1000):
    return [sequence[i:i+max_length] for i in range(0, len(sequence), max_length)]

# Function to find repeating amino acid sequences
def find_hetero_amino_acid_repeats(sequence):
    repeat_counts = defaultdict(int)

    # Iterate over all possible substring lengths
    for length in range(2, len(sequence) + 1):
        for i in range(len(sequence) - length + 1):
            substring = sequence[i:i+length]
            repeat_counts[substring] += 1

    # Filter out substrings that occur only once
    return {k: v for k, v in repeat_counts.items() if v > 1}

# Function to check and update repeats at boundaries
def check_boundary_repeats(fragments, final_repeats, overlap=50):
    """
    Check for repeating substrings that span across fragment boundaries
    and update the final repeats dictionary accordingly.
    
    Ensures that repeats are truly spanning both fragments.
    """
    for i in range(len(fragments) - 1):
        left_overlap = fragments[i][-overlap:] if len(fragments[i]) >= overlap else fragments[i]
        right_overlap = fragments[i + 1][:overlap] if len(fragments[i + 1]) >= overlap else fragments[i + 1]
        overlap_region = left_overlap + right_overlap  # Join both
        
        boundary_repeats = find_hetero_amino_acid_repeats(overlap_region)

        for substring, count in boundary_repeats.items():
            # Ensure substring spans across both fragments
            if any(aa in left_overlap for aa in substring) and any(aa in right_overlap for aa in substring):
                final_repeats[substring] += count  # Only add if spanning both fragments

    return final_repeats

# Function to find new repeats that only appear at fragmentation points
def find_new_boundary_repeats(fragments, final_repeats, overlap=50):
    """
    Identify new repeats that appear only at fragmentation boundaries
    and update the final dictionary accordingly.
    
    Ensures that:
    - The detected repeat spans across two fragments.
    - It is *not* already counted from individual fragments.
    """
    new_repeats = defaultdict(int)

    for i in range(len(fragments) - 1):
        left_overlap = fragments[i][-overlap:] if len(fragments[i]) >= overlap else fragments[i]
        right_overlap = fragments[i + 1][:overlap] if len(fragments[i + 1]) >= overlap else fragments[i + 1]
        overlap_region = left_overlap + right_overlap  # Join both
        
        boundary_repeats = find_hetero_amino_acid_repeats(overlap_region)

        for substring, count in boundary_repeats.items():
            # Ensure substring spans across both fragments
            if any(aa in left_overlap for aa in substring) and any(aa in right_overlap for aa in substring):
                # Only add if it's a *new* repeat (not already in final dictionary)
                if substring not in final_repeats:
                    new_repeats[substring] += count  

    return new_repeats

# Main function to process the protein sequence
def process_protein_sequence(sequence, overlap=50):
    fragments = fragment_protein_sequence(sequence)
    
    # Step 1: Find repeats in each fragment
    final_repeats = defaultdict(int)
    for fragment in fragments:
        fragment_repeats = find_hetero_amino_acid_repeats(fragment)
        for k, v in fragment_repeats.items():
            final_repeats[k] += v

    # Step 2: Check and update repeats at boundaries
    final_repeats = check_boundary_repeats(fragments, final_repeats, overlap)

    # Step 3: Find new repeats emerging at boundaries
    new_repeats = find_new_boundary_repeats(fragments, final_repeats, overlap)

    # Step 4: Merge new repeats into final dictionary
    for k, v in new_repeats.items():
        final_repeats[k] += v

    return final_repeats

# Example usage
if _name_ == "_main_":
    sequence_length = 12030  # Example length
    protein_sequence = generate_protein_sequence(sequence_length)
    
    final_repeats = process_protein_sequence(protein_sequence)
    
    print("Final Repeats Found:", final_repeats)