Title: Target Detection Tabulation

# Detection Tabulation

We make a tabulation of the detection process in which 
target stars or planets are found and a count of 
(hopefully) successful detections accumulates.
There are two tables, each restricted to a specific target family:

* All stars: Every star in the target list which was observed.

* Promoted stars: a list of stars that have been promoted into 
a category where they are eligible for spectral characterization.
The EXOSIMS scheduler controls promotion, and provides the
list at the end of the simulation in the `promoted_stars` variable.

All quantities shown are ensemble averages of detection-related
events across all DRM's in the ensemble - the ensemble mean.

All quantities after &plusmn; are the
[standard error of the mean](https://en.wikipedia.org/wiki/Standard_error).
That is, they are the population standard deviation divided by
the square root of ensemble size, so they will drop to zero
as ensemble size grows.
The ensemble size used for this scaling is shown in a note below the table.
The intent is to show how accurate the shown average is.

## Table Contents

In all cases, consideration is restricted to a population P of
stars (all observed stars, or promoted stars).

Every detection observation is an "Attempt," ending in either
Failure (`det_status` = 0, -1, or -2 in the simulation DRM),
or Success (`det_status = 1` in the DRM).

In general, "cumulative" means the total number of detections,
counting repeated detections of the same target. And
"unique" means that repeat detections, or failures, are not counted.

Specific notes on row-by-row properties:

- Attempts end in either failure or successful detection (`A = F + D`).
  <br>
  This carries through to cumulative numbers in the table, so that:
  <br>
 `Attempt_cume = Fails_cume + Detected_cume`
   <br>
   This property does not carry through to unique counts. 
   (E.g., two attempts on one target could end in one failure and one
   success, and all unique numbers would equal 1.)

- Failures also sub-divide, so that:
  <br>
  `Fails_Cume = Fails_IWA + Fails_OWA + Fails_SNR`
  <br>
	All of the reported subcategories are cumulative counts.

- Failure categories do not make sense for stars because there is no way to
  attribute a specific set of planet failures to the host star.
  They are blanked out in the table.

- `Detected(n)` for planets means:
  <br>
  "the number of planets having exactly `n`
  occurrences of (`det_status=1`) during the mission". 
  <br>
  There is also an overflow category, `Det(n > 4)`.

- `Detected(n)` for stars is less obvious. We choose it to mean:
  <br> 
  "the number of
  stars with any combination of planets having exactly `n` ocurrences of (`det_status=1`) during the
  mission". 
  <br>
  This is generally useful only for Earth-only SU's, or SU's that are
  restricted to one planet type.

- Because the `Detected(n)` counts partition the set of all detections, 
  we get this relationship:
  <br>
  `Detected_uniq = Det(n=1) + Det(n=2) + Det(n=3) + Det(n=4) + Det(n>4)`

- Note that the `Detected(n)` categories imply that `n` detections were cumulatively performed.
  So, provided the overflow class (`Detected(n>4)`) is empty, we also have:
  <br>
  `Detected_cume = 1*Det(n=1) + 2*Det(n=2) + 3*Det(n=3) + 4*Det(n=4)`






