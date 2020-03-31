# Changes in version 0.2.6 released on 31 Mar 2020
## New features
## Bug fixes
1. Fixed a bug that the number of generations for candidates would be traced to n+1 when tracegen=n. This bug is found by Mianyu Liu.

# Changes in version 0.2.5 released on 25 Feb 2020
## New features
## Bug fixes
1. The tidyped() does not work with trace='all' in [certain cases](https://github.com/luansheng/visPedigree/issues/2#issue-568599008)

# Changes in version 0.2.4.1 released on 24 Feb 2020
## New features
## Bug fixes
1. An unexpected column with the name as NA occured when a tidyped object is tidyed again using the tidyped()

# Changes in version 0.2.4 released on 12 June 2019
## New features
## Bug fixes
1. The data.table used as the input parameter 'ped' may be changed in tidyped() and visped().


# Changes in version 0.2.3 released on 05 Mar 2019
## New features
## Bug fixes
1. The generation number of individuals is not inferred rightly.

# Changes in version 0.2.2 released on 28 Jan 2019
## New features
## Bug fixes
1. The tidied pedigree will not include the candidates which are not in the Ind column of the origin pedigree when the cand parameter is not NULL.

# Changes in version 0.2.1 released on 17 Nov 2018
## New features
## Bug fixes
1. Repel the overlapping nodes due to very small differences (digits > 7) among x positions of nodes
