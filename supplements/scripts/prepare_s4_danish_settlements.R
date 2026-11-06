#!/usr/bin/env Rscript

# Supplement 4: derive the Danish Neolithic settlement dataset
#
# Input:  ../data/anlaeg_all.csv
#         Full export of all features ("anlæg") from the Danish national
#         register Fund og Fortidsminder (Slots- og Kulturstyrelsen),
#         downloaded 5 April 2018 (https://www.kulturarv.dk/fundogfortidsminder/).
#         The file is Latin-1 encoded.
#
# Output: ../data/anlaeg_settlements.csv (UTF-8)
#
# Selection:
#   1. datering == "Stenalder" (Stone Age)
#   2. feature type (anlaegstype) in the list of settlement-related types
#      below (settlement, houses, pits, hearths, cultural layers, wells,
#      fences, workshops, hunting/fishing sites, fortifications)
#   3. fra_aar >= -4200 and til_aar <= -1701, i.e. records dated within the
#      Neolithic (end dates in the source are inclusive; see Supplement 4)
#
# Note: shell middens (type 3188, "Skaldynge") and flint-knapping sites
# (15437, "Flinthuggeplads") are not included.
#
# Run from the supplements/ directory:
#
#   Rscript scripts/prepare_s4_danish_settlements.R
#
# By default the existing output file is not overwritten; the result is
# compared with it instead. Set OVERWRITE=true to replace it.
library(dplyr)

input_file <- "../data/anlaeg_all.csv"
output_file <- "../data/anlaeg_settlements.csv"
overwrite <- Sys.getenv("OVERWRITE", "false") == "true"

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

settlement_types <- c(
  3000,  # Bosættelse, uspec. undergruppe
  3028,  # Brønd
  3030,  # Bygning
  3081,  # Hegn/gærde
  3089,  # Hus (evt. med stald)
  3097,  # Ildsted
  3125,  # Kogegrube
  3131,  # Kulturlag
  3427,  # Affaldsgrube
  7030,  # Bygning
  7040,  # Dyreben
  7070,  # Grube (uspecificeret funktion)
  7097,  # Ildsted
  7113,  # Kildeanlæg
  7125,  # Kogegrube
  7131,  # Kulturlag
  7427,  # Affaldsgrube
  7515,  # Ovnanlæg (uspecificeret)
  7524,  # Stolpehul m. uvis funktion
  13000, # Primær fødeerhverv, uspec. undergruppe
  13131, # Kulturlag
  13433, # Fiskeplads
  13434, # Indsamling-/jagtplads
  15131, # Kulturlag
  15256, # Værksted
  15427, # Affaldsgrube
  26029  # Befæstning
)

raw <- readr::read_csv(
  input_file,
  locale = readr::locale(encoding = "latin1"),
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)

settlements <- raw |>
  mutate(anlaegstype = as.integer(anlaegstype), fra_aar = as.numeric(fra_aar),
         til_aar = as.numeric(til_aar)) |>
  filter(datering == "Stenalder", anlaegstype %in% settlement_types,
         fra_aar >= -4200, til_aar <= -1701)

message("Selected ", nrow(settlements), " of ", nrow(raw), " records.")

key_cols <- c("systemnr", "anlnr", "anlaegstype", "fra_aar", "til_aar")

if (file.exists(output_file) && !overwrite) {
  existing <- readr::read_csv(
    output_file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ) |>
    mutate(anlaegstype = as.integer(anlaegstype), fra_aar = as.numeric(fra_aar),
           til_aar = as.numeric(til_aar))

  new_keys <- settlements |> mutate(systemnr = as.character(as.integer(systemnr))) |> count(across(all_of(key_cols)))
  old_keys <- existing |> mutate(systemnr = as.character(as.integer(systemnr))) |> count(across(all_of(key_cols)))
  identical_selection <- nrow(settlements) == nrow(existing) &&
    nrow(anti_join(new_keys, old_keys, by = c(key_cols, "n"))) == 0

  message("Existing file ", output_file, " has ", nrow(existing), " records. ",
          if (identical_selection) "The selection is identical." else
          "The selection DIFFERS; set OVERWRITE=true to replace the file.")
} else {
  readr::write_csv(settlements, output_file)
  message("Written: ", output_file)
}
