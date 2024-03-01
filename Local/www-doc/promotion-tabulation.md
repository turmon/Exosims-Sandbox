Title: Target Promotion Tabulation

# Promotion Tabulation

We make a tabulation of a "promotion funnel" in which potential 
target stars or planets are found, promoted, and characterized.
There are two tables, each restricted to a specific target family:

* Deep-dive stars: a list of Hipparcos numbers for target stars of 
special interest, controlled by the `top_HIPs` script variable.
Typically this list is smaller than the list of all promoted stars below.

* Promoted stars: a list of stars that have been promoted into
a category where they are eligible for spectral characterization.
The EXOSIMS scheduler controls promotion, and provides the
list at the end of the simulation in the `promoted_stars` variable.

**It is helpful to visualize the table as a funnel**, in which stars are
promoted
to be eligible for spectral characterization (row 1),
characterizations are attempted (rows 2 and 3),
and some of them succeed (rows 4 and 5).

The same is true of planets: in order to be successfully characterized
(last rows) the star around the planet must be promoted (row 1).

## Table Contents

In all cases, consideration is restricted to a population P of
stars (deep-dive or promoted).
All quantities in the table are averages across the ensemble.

Below is an example table for an Earths-only scenario
as a visual reminder.

[comment]: # (markdown processor doesn't allow multiline table head) 
[comment]: # (pasted from a markdown table .../tbl/table-funnel.md)

| Status                   |                     Star |                  Planet |               Hab. Zone |                   Earth |
| ------------------------ | ------------------------ | ------------------------ | ------------------------ | ------------------------ |
| Promoted                 |      19.83 &plusmn; 0.38 |     21.94 &plusmn; 0.45 |     21.94 &plusmn; 0.45 |     21.94 &plusmn; 0.45 |
| Attempt: Cume            |      28.58 &plusmn; 1.23 |     31.85 &plusmn; 1.41 |     31.85 &plusmn; 1.41 |     31.85 &plusmn; 1.41 |
| Attempt: Unique          |      12.72 &plusmn; 0.32 |     14.46 &plusmn; 0.38 |     14.46 &plusmn; 0.38 |     14.46 &plusmn; 0.38 |
| Observed: Cume           |      28.58 &plusmn; 1.23 |     30.37 &plusmn; 1.28 |     30.37 &plusmn; 1.28 |     30.37 &plusmn; 1.28 |
| Observed: Unique         |      12.72 &plusmn; 0.32 |     13.98 &plusmn; 0.36 |     13.98 &plusmn; 0.36 |     13.98 &plusmn; 0.36 |


### Row 1: Promoted

Star
: How many stars were promoted. 
  If P is promoted stars, this is the ensemble average of `len(promoted_stars)`.
  If P is deep-dive stars, this will be the ensemble average of
  the intersection of `promoted_stars` and the deep-dive list.

Planet
:  How many planets orbit the promoted stars. See the note below.

Hab. Zone
: How many HZ planets orbit the promoted stars.

Earth
: How many Earthlike planets orbit the promoted stars.

Note: The planet counts in this row are counting
*how many planets are around promoted stars*.
They do not say anything about 
observing (detecting) any specific planets.
In the illustrated case, there was on average
about 1.1 planet per promoted star (21.94 planets, 19.83 stars).
This is consistent with an "Earths-only" planet population:
all promoted stars will have one planet, and some will have two 
or perhaps three.

This definition follows the "funnel" concept
mentioned above.
The promotion mechanism has reduced
the pool of stars eligible for characterization
from the full `TargetList` to just the promoted stars.
Correspondingly, the pool of possible planets is now reduced to
that in the right-hand three columns.

Depending on observation scheduling and success, the pool
of planets will be reduced further, as seen in the next rows.

### Rows 2-3: Attempted Spectral Characterizations 

First, the "Attempt: Cume" row.

Star
:   The ensemble average of the number of 
	attempted spectral characterizations of stars in P.
    An alternative wording: the average number of characterization
    visits to stars in P, counting revisits.

Planet
:   The total number of planets potentially observed by such visits and
    re-visits, regardless of success. If a star with 8 planets is visited 3
    times, this number will increase by 24.

Hab. Zone, Earth 
:   The total number of Habitable Zone (resp. Earthlike) planets
    potentially observed by such visits, regardless of success.

The "Attempt: Unique" row simply ignores revisits. It may never
exceed the above row.

Star
:   The ensemble average of stars in P with any attempted spectral
	characterization.
	Alternatively: the number of stars in P with *any* characterization
	visit at all.

Planet, Hab. Zone, Earth 
:   The total number of planets (resp. HZ or Earthlike planets)
	of stars in P potentially observed by such visits.
	In the case above where a star with 8 planets is visited 3 times,
	this number increases by only 8, because revisits are not counted.

Again, the "funnel" concept is evident:

* Our efficiency in scheduling the promoted stars is measured by the
	difference between these rows and the first row. This difference can be
	denominated in either Stars or Planets. 

* Observing opportunities spent on revisits are measured by the 
  difference between Cumulative versus Unique attempts.

In the scenario above, we visited about 2/3 of the promoted stars
(12.72 visits to 19.83 stars), and we typically revisited stars at
least once (28.58 total visits, 12.72 unique visits).

### Rows 4-5: Successful Spectral Characterizations

First, the "Observed: Cume" row.

In the descriptions below, "successful characterization" is
the generous interpretation: 
any partial or full spectral characterization, in any band.
 
Star
:   The count of all successful spectral characterizations
	of any planet of a star in P, counting repeat characterizations.
    An alternative wording is: the count of successful char visits 
	and re-visits to stars in P. 
	(Successful means that at least one planet around that
	star was characterized.)

Planet, Hab. Zone, Earth 
:   The count of all successful characterizations of planets 
	(resp., HZ or Earthlike planets) of stars in P.
	In this case, unlike for Attempts,
	the planet-by-planet success is taken in to account.
	So in the case above of the star with 8 planets, the particular
	characterization status of each planet will be added up over all repeat
	visits. 

Again, the "Observed: Unique" row ignores revisits, and it
may never exceed the above row.

Star
:   The count of stars in P with any successful spectral characterization.
	Alternatively: the count of stars in P with a successful
	characterization visit.
	
Planet, Hab. Zone, Earth 
:   The number of planets (resp. HZ or Earthlike planets)
	of stars in P 
	characterized at least once.
	
This segment of the "funnel" covers our ability to schedule
successful characterizations, e.g., at favorable phase angle
and satisfying keepout.




