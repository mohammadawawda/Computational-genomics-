# Library complexity

`library_complexity.py` reads the fragment files of sequencing libraries and counts how many times every fragment was seen. It fits a zero truncated Poisson to these counts, estimates how many different molecules the library really has, and checks the fit with a chi squared test. It also checks how the estimate changes with the fragment length and saves the plots.

`report.pdf` is my written report.

## How to run

Put the data in a folder named `data` next to the script, then run `python library_complexity.py`.
