# current-loss-analysis
Estimates current density loss and current density generated (both totals and bucketed) for solar cell given quantum efficiency data
#Author: Alec Mullen
#Company: Siva Power 

This code first partitions the QE curve at critical points using first and second derivatives based on rough estimates 
of where the curve should change given the chemical properties of each of the cell's layers. It then performs a best fit 
for each of these partitions. These best fits are convolved with the AM1.5G to find current generated and current loss.
Current generated and current loss are also given for each of the partitions, since they will give insight into the 
effects of absorption at each layer. The best fits used, electron current density curve vs AM1.5G, and bucketed current 
generated and current loss are all graphed in a JMP dashboard.

Run JSL Script on JMP Data Table with QE (quantum efficiency in %), WL (wavelength in nm), and Bandgap (in eV) columns.
Also uses Excel file w/ AM1.5G Photon flux data (WL in first column, Photons/ (nm * m^2 * s) in second) w/ lables in 
first row. You must change path names in JSL script first. It will talk to python to perform the analysis. Makes use of 
numpy, openpyxl, and matplotlib. matplotlib is not used when JMP is used, so it can be commented out or deleted. 

The JSL script talks to Python using the method Run Program() and communicates via stdin and stout sending and parsing
strings. To debug python errors, you must encapsulate << Read statements in a Show() statement to see tracebacks. It
might also be useful to write a JSL method for reading that performs the << Read and prints it out iff the string is 
an error message.

By changing global variable JMP to False in the python script, you can run the analysis using an Excel file instead.
So far, it will graph the best fit and print out total current loss only. Using matplotlib and the functions that are 
already written (specifically current_loss_all_data()) it should be trivial to plot or print out relevant data. You will
also have to change the path name for the Excel file used. the Excel file should have the same columns that the JMP
Table would need, but w/ lables in the first row.
