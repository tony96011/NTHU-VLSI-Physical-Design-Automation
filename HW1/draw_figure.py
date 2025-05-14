import pandas as pd
import matplotlib.pyplot as plt

def plot_data(csv_file, x_column, y_column, title):
    """Plots the effect of x_column on y_column and saves the figure."""
    # Load the data
    df = pd.read_csv(csv_file)

    # Create a figure
    plt.figure(figsize=(10, 6))

    # Plot the data
    plt.plot(df[x_column], df[y_column], label=y_column, marker="o", linestyle="-")

    # Labels and legend
    plt.xlabel(x_column)
    plt.ylabel(y_column)
    plt.title(title)
    plt.legend()
    plt.grid()

    # Generate output filename dynamically
    output_file = f"{x_column}_{y_column.lower()}.png"
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Figure saved as {output_file}")

# Define datasets and plots
cu_csv = "results_cu.csv"
cp_csv = "results_cp.csv"

# Plot CU-based figures
plot_data(cu_csv, "CU", "Total_Chip_Area", "Effect of CU on Total Chip Area")
plot_data(cu_csv, "CU", "Total_Wire_Length", "Effect of CU on Total Wire Length")
plot_data(cu_csv, "CU", "Slack_Time", "Effect of CU on Slack Time")

# Plot CP-based figures
plot_data(cp_csv, "CP", "Total_Chip_Area", "Effect of CP on Total Chip Area")
plot_data(cp_csv, "CP", "Total_Wire_Length", "Effect of CP on Total Wire Length")
plot_data(cp_csv, "CP", "Slack_Time", "Effect of CP on Slack Time")
