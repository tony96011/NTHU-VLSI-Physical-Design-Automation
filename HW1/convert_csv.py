import csv

# Define file names
txt_file = "results_summary.txt"  # Replace with your actual TXT file name
csv_file = "results.csv"  # Output CSV file

# Read TXT and write to CSV
with open(txt_file, "r", encoding="utf-8") as infile, open(csv_file, "w", newline="", encoding="utf-8") as outfile:
    lines = infile.readlines()

    # Extract header
    header = lines[0].split()  # Splitting on whitespace
    separator_index = 1  # Line index where separator (----) is found

    # Write to CSV
    writer = csv.writer(outfile)
    writer.writerow(header)  # Write header

    # Write data rows
    for line in lines[separator_index + 1:]:
        row = line.split()
        writer.writerow(row)

print(f"Conversion complete! Saved as {csv_file}")
