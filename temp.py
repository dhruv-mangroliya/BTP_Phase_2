import random
import xlsxwriter
from io import BytesIO
from collections import defaultdict
import streamlit as st
import pandas as pd

# Function to generate a random protein sequence of given length
def generate_protein_sequence(length):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # 20 standard amino acids
    return ''.join(random.choices(amino_acids, k=length))

# Function to fragment the protein sequence into chunks of max length 1000
def fragment_protein_sequence(sequence, max_length=1000):
    return [sequence[i:i+max_length] for i in range(0, len(sequence), max_length)]

# Function to find heterorepeats in the protein sequence
def find_hetero_amino_acid_repeats(sequence):
    n = len(sequence)
    freq = defaultdict(int)

    # Iterate through all possible lengths of substrings (2 to N)
    for length in range(2, n + 1):
        # Use a sliding window to find substrings of this length
        for i in range(n - length + 1):
            substring = sequence[i:i + length]

            # Check if all amino acids in the substring are different
            if len(set(substring)) == len(substring):  # All amino acids are unique
                freq[substring] += 1

    # Filter out repeats with frequency 1 or length 1
    return {repeat: count for repeat, count in freq.items() if len(repeat) > 1 and count > 1}

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

# Function to process the protein sequence and return final heterorepeats
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

    return final_repeats

# Function to process a single Excel sheet and return its analysis
def process_excel(excel_data):
    sequence_data = []
    all_heterorepeats = defaultdict(int)  # Track all heterorepeats and their counts

    for sheet_name in excel_data.sheet_names:
        df = excel_data.parse(sheet_name)
        if len(df.columns) < 3:
            st.error(f"Error: The sheet '{sheet_name}' must have at least three columns: ID, Protein Name, Sequence")
            return None, None

        for _, row in df.iterrows():
            entry_id = str(row[0])
            protein_name = str(row[1])
            sequence = str(row[2]).replace('"', '').replace(' ', '')
            freq = find_hetero_amino_acid_repeats(sequence)  # Returns a dictionary with repeats and their counts
            sequence_data.append((entry_id, protein_name, freq))

            # Update the main heterorepeats dictionary with counts
            for repeat, count in freq.items():
                all_heterorepeats[repeat] += count

    return all_heterorepeats, sequence_data

# Function to generate and download Excel workbook with separate sheets for each input file
def create_excel(sequences_data, heterorepeats, filenames):
    output = BytesIO()
    workbook = xlsxwriter.Workbook(output, {'in_memory': True})

    # Filter out repeats with frequency 1 or length 1
    valid_heterorepeats = {repeat: count for repeat, count in heterorepeats.items() if len(repeat) > 1 and count > 1}

    # Iterate through sequences data grouped by filenames and create separate sheets
    for file_index, file_data in enumerate(sequences_data):
        filename = filenames[file_index]
        worksheet = workbook.add_worksheet(filename[:31])  # Limit sheet name to 31 characters

        # Write the header for the current file
        worksheet.write(0, 0, "Entry ID")
        worksheet.write(0, 1, "Protein Name")
        col = 2
        for repeat in sorted(valid_heterorepeats):
            worksheet.write(0, col, repeat)
            col += 1

        # Write data for each sequence in the current file
        row = 1
        for entry_id, protein_name, freq in file_data:
            worksheet.write(row, 0, entry_id)
            worksheet.write(row, 1, protein_name)
            col = 2
            for repeat in sorted(valid_heterorepeats):
                worksheet.write(row, col, freq.get(repeat, 0))
                col += 1
            row += 1

    workbook.close()
    output.seek(0)
    return output

# Streamlit UI components
st.title("Protein Heterorepeat Analysis")

# Step 1: Upload Excel Files
uploaded_files = st.file_uploader("Upload Excel files", accept_multiple_files=True, type=["xlsx"])

# Step 2: Process files and display results
if uploaded_files:
    all_heterorepeats = defaultdict(int)
    all_sequences_data = []
    filenames = []

    for file in uploaded_files:
        excel_data = pd.ExcelFile(file)
        heterorepeats, sequence_data = process_excel(excel_data)
        if heterorepeats is not None:
            all_heterorepeats.update(heterorepeats)
            all_sequences_data.append(sequence_data)
            filenames.append(file.name)

    if all_sequences_data:
        st.success(f"Processed {len(uploaded_files)} files successfully!")

        # Step 3: Generate and download the Excel report
        excel_file = create_excel(all_sequences_data, all_heterorepeats, filenames)

        # Download the Excel file
        st.download_button(
            label="Download Excel file",
            data=excel_file,
            file_name="protein_heterorepeat_results.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

        # Step 4: Display summary table
        if st.checkbox("Show Results Table"):
            # Convert the sequences data into a DataFrame for easy display
            rows = []
            for file_index, file_data in enumerate(all_sequences_data):
                filename = filenames[file_index]
                for entry_id, protein_name, freq in file_data:
                    row = {"Filename": filename, "Entry ID": entry_id, "Protein Name": protein_name}
                    row.update({repeat: freq.get(repeat, 0) for repeat in sorted(all_heterorepeats)})
                    rows.append(row)

            result_df = pd.DataFrame(rows)
            st.dataframe(result_df)
