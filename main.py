# Flow Cell Performance Dashboard Built from Public ONT Run Telemetry
# Krutik Peddakotla -- https://www.linkedin.com/in/krutik-peddakotla/

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

file = pd.read_csv("sequencing_summary_PAG51911_6e676ba2.txt", sep="\t")
df = pd.DataFrame(file)
df = df.sort_values(by=["start_time"], ascending=True)
df.drop(columns=["num_events","template_start","template_duration",
                 "filename_fastq", "filename_fast5"], inplace=True)

# Validation #
df["end_time"] = ((df["start_time"] + df["duration"]) / 3600)
df["speed"] = (df["sequence_length_template"] / df["duration"])
print(df.head().to_string())

print(f"Checking for unique run ids: {df["run_id"].value_counts()} \n")
print(f"Max Channel #: {df["channel"].max()}")
print(f"Num Unique Channel #s: {df["channel"].nunique()} \n")
print(f"Min Start Time: {df["start_time"].min()}")
print(f"Max Start Time: {df["start_time"].max()}")
print(f"Start Time Range: {df["start_time"].max() - df["start_time"].min()} \n")
print(f"Checking for unique mux values: {df["mux"].value_counts()} \n")
print(f"Checking for end reasons: {df["end_reason"].value_counts()} \n")
print(f"Checking for number of passes: {df["passes_filtering"].value_counts()} \n")


## Distributions & Graphing ##
print(f"Read Length Distributions: \n{df["sequence_length_template"].describe()} \n")
print(f"Q-Score Distributions: \n{df["mean_qscore_template"].describe()} \n")
print(f"Duration Distributions: \n{df["duration"].describe()} \n")

sns.set_theme(style="whitegrid", context="notebook")
plt.figure(figsize=(10, 6))
plt.subplot(3, 2, 1)

readlength_plot = sns.histplot(data=df, x="sequence_length_template", stat="density", common_norm=False,
                         bins=100, edgecolor=None, alpha=0.75, log_scale=True, element="step",
                         hue="passes_filtering", multiple="layer")
plt.xlabel("Sequence Length", fontsize=12)
plt.ylabel("Density (per log10 unit)", fontsize=12)
plt.title("Sequence Length Distributions (Normalized)", fontsize=14, pad=12)

plt.subplot(3, 2, 2)
qscore_plot = sns.histplot(data=df, x="mean_qscore_template", stat="density", common_norm=False,
                         bins=100, edgecolor=None, alpha=0.75, element="step",
                         hue="passes_filtering", multiple="layer")
plt.xlabel("Q-Score", fontsize=12)
plt.ylabel("Density", fontsize=12)
plt.title("Q-Score Distribution (Normalized)", fontsize=14, pad=12)

plt.subplot(3, 2, 3)
duration_plot = sns.histplot(data=df, x="duration", stat="density", common_norm=False,
                         bins=100, edgecolor=None, alpha=0.75, log_scale=True, element="step",
                         hue="passes_filtering", multiple="layer")
plt.xlabel("Duration", fontsize=12)
plt.ylabel("Density (per log10 unit)", fontsize=12)
plt.title("Duration Distribution (Normalized)", fontsize=14, pad=12)

plt.subplot(3, 2, 4)
hues=["mux_change","unblock_mux_change"]
filtered_df = df[df["end_reason"].isin(hues)]
endtime_plot = sns.histplot(data=filtered_df, x="end_time",
                         bins=100, edgecolor=None, alpha=0.75, element="step",
                         hue="end_reason", multiple="layer")
plt.yscale("log")
plt.xlabel("End Time (hours)", fontsize=12)
plt.ylabel("Count (per log10 unit)", fontsize=12)
plt.title("End Time", fontsize=14, pad=12)

plt.subplot(3, 2, 5)
speed_distribution = sns.histplot(data=df, x="speed",
                         bins=100, edgecolor=None, alpha=0.75)
plt.xlabel("Speed", fontsize=12)
plt.ylabel("Count", fontsize=12)
plt.title("Speed Distribution", fontsize=14, pad=12)
plt.tight_layout()

## Determine minimum q-score threshold for pass/fail ##
#print(df[df["passes_filtering"] == True].sort_values(by=["mean_qscore_template"], ascending=True).head().to_string())
#print(df[df["passes_filtering"] == False].sort_values(by=["mean_qscore_template"], ascending=False).head().to_string())

df.to_parquet("summary.parquet.gzip",index=False, compression="gzip")
#plt.show()