import os
import streamlit as st
import pandas as pd
import xlsxwriter
from io import BytesIO
from collections import defaultdict

# Function to find heterorepeats in the protein sequence
def find_heterorepeats(protein):
    n = len(protein)
    freq = defaultdict(int)

    # Iterate through all possible lengths of substrings (2 to N)
    for length in range(2, n + 1):
        # Use a sliding window to find substrings of this length
        for i in range(n - length + 1):
            substring = protein[i:i + length]

            # Check if all amino acids in the substring are different
            if len(set(substring)) == len(substring):  # All amino acids are unique
                freq[substring] += 1

    return freq

# Function to process a single Excel sheet and return its analysis
def process_excel(excel_data):
    heterorepeats = set()
    sequence_data = []

    for sheet_name in excel_data.sheet_names:
        df = excel_data.parse(sheet_name)
        if len(df.columns) < 3:
            st.error(f"Error: The sheet '{sheet_name}' must have at least three columns: ID, Protein Name, Sequence")
            return None, None

        for _, row in df.iterrows():
            entry_id = str(row[0])
            protein_name = str(row[1])
            sequence = str(row[2]).replace('"', '').replace(' ', '')
            freq = find_heterorepeats(sequence)
            sequence_data.append((entry_id, protein_name, freq))
            heterorepeats.update(freq.keys())  # Collect unique heterorepeats

    return heterorepeats, sequence_data

# Function to generate and download Excel workbook with separate sheets for each input file
def create_excel(sequences_data, heterorepeats, filenames):
    output = BytesIO()
    workbook = xlsxwriter.Workbook(output, {'in_memory': True})

    # Iterate through sequences data grouped by filenames and create separate sheets
    for file_index, file_data in enumerate(sequences_data):
        filename = filenames[file_index]
        worksheet = workbook.add_worksheet(filename[:31])  # Limit sheet name to 31 characters

        # Write the header for the current file
        worksheet.write(0, 0, "Entry ID")
        worksheet.write(0, 1, "Protein Name")
        col = 2
        for repeat in sorted(heterorepeats):
            worksheet.write(0, col, repeat)
            col += 1

        # Write data for each sequence in the current file
        row = 1
        for entry_id, protein_name, freq in file_data:
            worksheet.write(row, 0, entry_id)
            worksheet.write(row, 1, protein_name)
            col = 2
            for repeat in sorted(heterorepeats):
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
    all_heterorepeats = set()
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
