[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/dunn&file=dunn.m)

# dunn

Stepdown Dunn procedure for nonparametric multiple comparisons between independent groups, implemented in MATLAB.

This function performs Dunn's stepdown procedure using rank-based statistics to compare multiple independent samples, with an optional control group and structured outputs.

## MATLAB Online

You can open and run the function directly in MATLAB Online via:

https://matlab.mathworks.com/open/github/v1?repo=dnafinder/dunn&file=dunn.m

## Repository

GitHub repository:

https://github.com/dnafinder/dunn

Main file:

- dunn.m – implementation of Dunn's stepdown procedure, with Name–Value options and structured outputs.

## Syntax

The function uses a single intuitive calling style: one data vector per group, followed by optional Name–Value pairs.

Basic forms:

- dunn(X1, X2, ..., XK)
- dunn(X1, X2, ..., XK, 'Ctrl', CTRL, 'Display', DISPLAY)
- stats = dunn(X1, X2, ..., XK, ...)
- [stats, results] = dunn(X1, X2, ..., XK, ...)

where:

- X1, ..., XK are numeric vectors (row or column) containing the observations for each of the K groups.
- CTRL is a logical-like flag indicating whether the first group X1 is to be treated as a control group.
- DISPLAY is a logical-like flag controlling command-window output.

Logical-like values accepted for CTRL and DISPLAY include:

- logical: true, false
- numeric: 1, 0
- char/string (case-insensitive): "on", "off", "yes", "no", "true", "false"

## Inputs

### Required inputs

- X1, ..., XK  
  - Numeric vectors (row or column), real, finite, non-NaN, non-empty.  
  - Each Xi contains the data for the i-th group.  
  - Group sizes may be unequal.  
  - At least two groups are required (K ≥ 2).

### Name–Value pair arguments

- 'Ctrl'  
  - Logical-like flag.  
  - If true, the first group X1 is treated as a control group and only comparisons versus this control group are performed (K − 1 comparisons).  
  - If false (default), all pairwise comparisons are performed with stepdown adjustment.

- 'Display'  
  - Logical-like flag.  
  - If true (default), the function prints a summary of groups and pairwise comparisons to the command window.  
  - If false, the function runs silently and only returns outputs.

## Outputs

The function supports zero, one, or two output arguments.

### No output

dunn(X1, X2, ..., XK)  
dunn(X1, X2, ..., XK, 'Ctrl', true, 'Display', false)

If no output arguments are requested:

- When Display is true (default), the function prints:
  - a group summary table,
  - the ties factor,
  - the pairwise comparison results table.
- When Display is false, no output is printed and nothing is returned.

### One output

stats = dunn(X1, X2, ..., XK, ...)

Returns a structure stats containing group-level summary statistics and global procedure parameters.

### Two outputs

[stats, results] = dunn(X1, X2, ..., XK, ...)

Returns both:

- stats – summary structure
- results – table of pairwise comparisons

### STATS structure

Fields:

- groups      – K×1 vector of group indices (1, 2, ..., K), in the order of input arguments.
- N           – K×1 vector of sample sizes per group.
- sumRanks    – K×1 vector of sums of ranks for each group.
- meanRanks   – K×1 vector of mean ranks for each group.
- tiesFactor  – scalar equal to 2 * tieadj, where tieadj is the tie adjustment returned by tiedrank.
- varFactor   – scalar f, the variance factor used in the denominator of the Q statistics, adjusted for ties.
- totalN      – total number of observations across all groups.
- alphaGlobal – global significance level (fixed at 0.05 in this implementation).
- alphaEff    – effective significance level used after Sidak adjustment for multiple comparisons.
- critQ       – critical Q value corresponding to alphaEff under the normal approximation.
- k           – number of groups.
- ctrl        – logical flag indicating whether a control group was used.

### RESULTS table

A table with one row per (potential) pairwise comparison. Columns:

- Comparison – text label indicating the pair of groups (for example "3-1").
- Q_value    – Dunn's Q statistic for that pair.
- Crit_Q     – critical Q value used for the procedure (same across rows).
- Comment    – textual decision on the null hypothesis, such as:
  - "Reject H0"
  - "Fail to reject H0"
  - "No comparison made" (for comparisons skipped by stepdown logic; H0 is accepted by stepdown rationale).

## Method

All observations from all groups are pooled and ranked using MATLAB's tiedrank function, with standard correction for ties.

For each group:

- The sum of ranks and mean rank are computed.
- Let:
  - Ntot be the total number of observations across all groups.
  - tieadj be the tie adjustment returned by tiedrank.

The variance factor f used in the denominator of the Q statistic is:

f = Ntot (Ntot + 1) / 12 − tieadj / [6 (Ntot − 1)]

For each pair of groups i and j, Dunn's Q statistic is:

Q = |mean_rank_i − mean_rank_j| / sqrt( f * (1/N_i + 1/N_j) )

where N_i and N_j are the group sample sizes.

### Multiple-comparison control and stepdown

Global significance level is set to 0.05.

#### No control group ('Ctrl', false)

- All K(K − 1)/2 pairwise comparisons are considered.
- A Sidak adjustment is applied so that the familywise confidence level is 0.95:
  - alphaEff = 1 − (1 − 0.05)^(1 / kstar)
  - kstar = K(K − 1)/2
- The critical Q value is obtained from the standard normal distribution via erfcinv.
- Groups are ordered by descending mean rank.
- Stepdown logic:
  - Comparisons proceed from the most extreme differences inward.
  - Once a non-significant comparison is encountered at a given step, intermediate comparisons become unnecessary and are labeled "No comparison made", with H0 accepted by stepdown reasoning.

#### With control group ('Ctrl', true)

- Only comparisons between the control (first group X1) and each of the remaining groups are performed (K − 1 comparisons).
- Sidak adjustment is applied for K − 1 comparisons:
  - alphaEff = 1 − (1 − 0.05)^(1 / (K − 1))
- The control group remains the first group; no sorting by mean rank is applied in this branch.

## Examples

### Example 1 – No control group, printed output

x1 = [7.68 7.69 7.70 7.70 7.72 7.73 7.73 7.76];  
x2 = [7.71 7.73 7.74 7.74 7.78 7.78 7.80 7.81];  
x3 = [7.74 7.75 7.77 7.78 7.80 7.81 7.84];  
x4 = [7.71 7.71 7.74 7.79 7.81 7.85 7.87 7.91];

dunn(x1, x2, x3, x4);

The function prints group-level summary statistics, ties factor, and the stepdown pairwise comparison table, including Q values and critical value.

### Example 2 – First group as control, silent run with outputs

[stats, results] = dunn(x1, x2, x3, x4, 'Ctrl', true, 'Display', false);

No output is printed. stats contains the group summary and procedure parameters, and results contains the pairwise comparison statistics and decisions versus the control.

## Requirements

- MATLAB R14 or later.
- The tiedrank function is required (Statistics and Machine Learning Toolbox in older MATLAB releases; integrated in many newer ones).
- No additional toolboxes are strictly required beyond those providing tiedrank and erfcinv.

## Usage Notes

- The legacy syntax using a pooled data vector X and a grouping vector G is not supported in this version. Users are expected to provide one data vector per group.
- Groups are internally labeled as 1, 2, ..., K in the order of the input arguments.
- The function is designed for independent groups; it is not intended for repeated-measures or paired designs.

## Citation

If you use this function in a scientific publication, please cite it as:

Cardillo G. (2006). Dunn's Test: a procedure for multiple nonparametric comparisons. Available on GitHub: https://github.com/dnafinder/dunn

## License

The code is provided as-is, without any explicit warranty. Please refer to the repository for licensing details if a LICENSE file is present.
