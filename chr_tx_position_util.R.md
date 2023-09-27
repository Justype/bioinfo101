# chr_tx_position_util.R

This is a script to convert chromosome position to transcript position.

You will notice there is a parameter called `is_5prime`. This is because R indexes are different from other programming languages. For the ease of slicing, the 5' indices of annotation are 1 nt less than the actual, like `exonStarts`, `cdsStarts`. However, the position in SAM/BAM file is the actual. So `is_5prime` aims to fix this problem.

| Function                                  | Description                                                |
| :---------------------------------------- | :--------------------------------------------------------- |
| [isInFrame](#isinframe)                   | check if the position is in the exon                       |
| [getTxPositions](#gettxpositions)         | get transcript positions based on chromosome positions     |
| [getFragTxStartEnd](#getfragtxstartend)   | get start & end indices based on annotation                |
| [getChrPositions](#getchrpositions)       | get chromosome positions based on transcript positions     |
| [getFragChrStartEnd](#getfragchrstartend) | get start & end indices based on annotation and transcript |

It is worth noting that the 5' end of the chromosome position is 1 nt less than the actual value, while it is the actual on the transcript.

# Functions

## isInFrame

```r
isInFrame = function(exon_starts, exon_ends, chr_position, is_5prime = F)
```

| Parameters     | Description                                    |
| :------------- | :--------------------------------------------- |
| `exon_starts`  | a vector of `exonStarts`                       |
| `exon_ends`    | a vector of `exonEnds`                         |
| `chr_position` | a chromosome position you want to check        |
| `is_5prime`    | is at 5' end of **chromosome** (e.g. cdsStart) |

| Returns | Condition              |
| :------ | :--------------------- |
| `TRUE`  | It is in one exon      |
| `FALSE` | It is not in all exons |

## getTxPositions

```r
getTxPositions = function(strand, exon_starts, exon_ends, chr_positions, is_5prime = F, is_check_frame = F)
```

| Parameters       | Description                                          |
| :--------------- | :--------------------------------------------------- |
| `strand`         | `+` or `-` strand                                    |
| `exon_starts`    | a vector of `exonStarts`                             |
| `exon_ends`      | a vector of `exonEnds`                               |
| `chr_positions`  | a vector of chromosome positions you want to convert |
| `is_5prime`      | is at 5' end of chromosome (e.g. cdsStart)           |
| `is_check_frame` | `TRUE` if is not annotation data                     |

Return: converted positions on transcript

## getFragTxStartEnd

```r
getFragTxStartEnd = function(strand, exon_starts, exon_ends, chr_start, chr_end, is_sam = F)
```

| Parameters    | Description                                   |
| :------------ | :-------------------------------------------- |
| `strand`      | `+` or `-` strand                             |
| `exon_starts` | a vector of `exonStarts`                      |
| `exon_ends`   | a vector of `exonEnds`                        |
| `chr_start`   | fragment 5' start (1 nt less if is_sam == F)  |
| `chr_end`     | fragment 3' end index (actual)                |
| `is_sam`      | `TRUE` if position is from SAM/BAM file <br> `FALSE` if positions from annotation like UCSC table |

Return: `c(tx_start_index, tx_end_index)` (Both index are the actual)

## getChrPositions

```r
getChrPositions = function(strand, exon_starts, exon_ends, tx_positions, is_5prime = F)
```

| Parameters     | Description                                          |
| :------------- | :--------------------------------------------------- |
| `strand`       | `+` or `-` strand                                    |
| `exon_starts`  | a vector of `exonStarts`                             |
| `exon_ends`    | a vector of `exonEnds`                               |
| `tx_positions` | a vector of transcript positions you want to convert |
| `is_5prime`    | is at 5' end of **transcript**                       |

Return: converted positions on chromosome (1 nt less if on 5' side of the chromosome)

## getFragChrStartEnd

```r
getFragChrStartEnd = function(strand, exon_starts, exon_ends, tx_start, tx_end)
```

| Parameters    | Description                                  |
| :------------ | :------------------------------------------- |
| `strand`      | `+` or `-` strand                            |
| `exon_starts` | a vector of `exonStarts`                     |
| `exon_ends`   | a vector of `exonEnds`                       |
| `tx_start`    | fragment 5' start (1 nt less if is_sam == F) |
| `tx_end`      | fragment 3' end index (actual)               |

Return: `c(chr_start_index, chr_end_index)` (1 nt less if on 5' side of the chromosome)
