library(dplyr)
library(ggplot2)
library(lubridate)
library(riem)
library(gridExtra)
library(cowplot)
library(tidyr)

setwd('/Users/precious/downloads/ASOS-PAPER')
asos_utc=read.csv('asos_utc.csv')

asos_utc=read.csv('test-asos.csv')

# simplest: alphabetical by station
asos_utc_sorted <- asos_utc %>% arrange(station)

library(dplyr)

asos_data_no_sa <- asos_utc_sorted %>%
  filter(!toupper(station) %in% "SA") %>%
  droplevels()

unique(asos_data_no_sa$station)

#------------------------------------calculate RH and add to the data------------------------------
library(dplyr)
library(lubridate)
library(stringr)

# --- helpers ---
safe_num <- function(x) suppressWarnings(as.numeric(x))

decode_M2 <- function(v){
  out <- rep(NA_real_, length(v))
  neg <- str_detect(v, "^M\\d{1,2}$")
  pos <- str_detect(v, "^\\d{1,2}$")
  out[neg] <- -safe_num(str_sub(v[neg], 2))
  out[pos] <-  safe_num(v[pos])
  out
}

# Parse T/Td: prefer "dd/dd" or "Mdd/Mdd"; fallback to RMK TsnTTTsdTdTdTd (tenths °C)
extract_T_Td <- function(metar_chr){
  n <- length(metar_chr)
  T_met  <- rep(NA_real_, n)
  Td_met <- rep(NA_real_, n)
  
  # 1) dd/dd with word boundaries (won’t match "1/2SM")
  m1 <- str_match(metar_chr, "(?<=\\s|^)(M?\\d{1,2})/(M?\\d{1,2})(?=\\s|$)")
  ok1 <- !is.na(m1[,1])
  if (any(ok1)) {
    T_met[ok1]  <- decode_M2(m1[ok1,2])
    Td_met[ok1] <- decode_M2(m1[ok1,3])
  }
  
  # 2) RMK T00000000 (tenths °C)
  m2 <- str_match(metar_chr, "\\bT(\\d{8})\\b")
  idx2 <- which(!is.na(m2[,2]))
  if (length(idx2)){
    v   <- m2[idx2,2]
    sT  <- as.integer(substr(v,1,1)); TTT <- as.integer(substr(v,2,4))
    sD  <- as.integer(substr(v,5,5)); DDD <- as.integer(substr(v,6,8))
    Tval <- (TTT/10) * ifelse(sT==1, -1, 1)
    Dval <- (DDD/10) * ifelse(sD==1, -1, 1)
    
    fill_idx <- idx2[is.na(T_met[idx2])]
    if (length(fill_idx)){
      sel <- is.na(T_met[idx2])
      T_met[fill_idx]  <- Tval[sel]
      Td_met[fill_idx] <- Dval[sel]
    }
  }
  tibble(T_met, Td_met)
}

# Magnus (Alduchov–Eskridge constants)
mag_rh <- function(Tc, Tdc){
  # RH = 100 * es(Td)/es(T)
  100 * exp((17.625*Tdc)/(243.04+Tdc) - (17.625*Tc)/(243.04+Tc))
}

# --- pipeline ---
parsed <- extract_T_Td(asos_data_no_sa$metar)

asos_data8 <- asos_data_no_sa %>%
  mutate(
    valid          = as.POSIXct(valid, tz = "UTC"),
    local_time     = with_tz(valid, tzone = "America/Los_Angeles"),
    date           = as.Date(local_time),
    year           = year(local_time),
    month          = month(local_time, label = TRUE, abbr = TRUE),
    hour           = hour(local_time),
    wind_speed_mps = safe_num(sknt) * 0.514444,
    visibility_km  = safe_num(vsby) * 1.60934
  ) %>%
  bind_cols(parsed) %>%
  mutate(
    # prefer METAR-parsed T/Td; fallback to existing tmpc/dwpc; else NA
    T_use  = coalesce(T_met,  safe_num(tmpc)),
    Td_use = coalesce(Td_met, safe_num(dwpc)),
    RH_calc = if_else(is.na(T_use) | is.na(Td_use),
                      NA_real_,
                      round(mag_rh(T_use, Td_use), 2)),
    rh_source = case_when(
      !is.na(T_met) & !is.na(Td_met) ~ "metar_pair_or_Tgroup",
      is.na(T_met) & !is.na(tmpc) & !is.na(dwpc) ~ "tmpc_dwpc_fallback",
      TRUE ~ "NA"
    ),
    month = factor(month, levels = month.abb, ordered = TRUE)
  ) %>%
  select(-T_met, -Td_met, -T_use, -Td_use)


asos_data9 <- asos_data8 %>%
  filter(year >= 2005, year <= 2024)

asos_data19 <- asos_data8 %>%
  filter(year >= 2005, year <= 2024)

dust_events_CV7 <- asos_data19 %>%
  filter(
    # Dust codes in either field
    grepl("\\b(DU|BLDU|DS|SS|VCBLDU)\\b", wxcodes, ignore.case = TRUE) |
      grepl("\\b(DU|BLDU|DS|SS|VCBLDU)\\b", metar,   ignore.case = TRUE) |
      
      # HZ (no FU) with thresholds, from wxcodes
      (wind_speed_mps >= 6 & visibility_km <= 10 & RH_calc < 70 &
         grepl("\\bHZ\\b", wxcodes, ignore.case = TRUE) &
         !grepl("\\bFU\\b", wxcodes, ignore.case = TRUE)) |
      
      # HZ (no FU) with thresholds, from metar
      (wind_speed_mps >= 6 & visibility_km <= 10 & RH_calc < 70 &
         grepl("\\bHZ\\b", metar,   ignore.case = TRUE) &
         !grepl("\\bFU\\b", metar,   ignore.case = TRUE))
  )


dust_events_cv7 <- dplyr::filter(dust_events_CV7, is.na(wind_speed_mps) | wind_speed_mps > 0)



dust_events_CV <- asos_data9 %>%
  filter(
    # Dust codes in either field
    grepl("\\b(DU|BLDU|DS|SS|VCBLDU)\\b", wxcodes, ignore.case = TRUE) |
      grepl("\\b(DU|BLDU|DS|SS|VCBLDU)\\b", metar,   ignore.case = TRUE) |
      
      # HZ (no FU) with thresholds, from wxcodes
      (wind_speed_mps >= 6 & visibility_km <= 10 & RH_calc < 70 &
         grepl("\\bHZ\\b", wxcodes, ignore.case = TRUE) &
         !grepl("\\bFU\\b", wxcodes, ignore.case = TRUE)) |
      
      # HZ (no FU) with thresholds, from metar
      (wind_speed_mps >= 6 & visibility_km <= 10 & RH_calc < 70 &
         grepl("\\bHZ\\b", metar,   ignore.case = TRUE) &
         !grepl("\\bFU\\b", metar,   ignore.case = TRUE))
  )


dust_events_cv5 <- dplyr::filter(dust_events_CV, is.na(wind_speed_mps) | wind_speed_mps > 0)

#---------CHECK OCTOBER 30-31 2019 [BFL]--CHECK THIS EVENT



# --------------------------------------------
# Daily PWC "ONLY" logic + ANY/MULTIPLE counts
# Input: dust_events5 with columns:
#   station, local_time (POSIXct, local tz), wxcodes, metar
# Output: one row per (station, date)
#   - *_only = count of reports where exactly that code appears (and no other of the set)
#   - ANY_PWC = count of reports with ≥1 of the codes
#   - MULTIPLE_PWC = count of reports with ≥2 of the codes
#   - *_flag, ANY_PWC_flag, MULTIPLE_PWC_flag = 0/1 day-level presence
# --------------------------------------------

library(dplyr)
library(stringr)
library(lubridate)

daily_pwc5 <- dust_events_cv5 %>%
  mutate(
    date   = as.Date(local_time),
    # combine both fields so we can search once
    txt_wx = paste(coalesce(as.character(wxcodes), ""),
                   coalesce(as.character(metar),   "")),
    
    # Base detections (strict tokens, case-insensitive)
    has_DU     = str_detect(txt_wx, regex("\\bDU\\b",     ignore_case = TRUE)),
    has_BLDU   = str_detect(txt_wx, regex("\\bBLDU\\b",   ignore_case = TRUE)),
    has_VCBLDU = str_detect(txt_wx, regex("\\bVCBLDU\\b", ignore_case = TRUE)),
    has_DS     = str_detect(txt_wx, regex("\\bDS\\b",     ignore_case = TRUE)),
    has_SS     = str_detect(txt_wx, regex("\\bSS\\b",     ignore_case = TRUE)),
    has_HZ     = str_detect(txt_wx, regex("\\bHZ\\b",     ignore_case = TRUE))
  ) %>%
  # Ensure VCBLDU does not double-count as BLDU
  mutate(
    has_BLDU_plain = has_BLDU & !has_VCBLDU,
    
    # Count how many of the codes appear in the SAME report
    n_pwc = rowSums(
      cbind(has_DU, has_BLDU_plain, has_DS, has_SS, has_VCBLDU, has_HZ),
      na.rm = TRUE
    ),
    
    # ONLY = exactly one code present in that report
    DU_only      = has_DU     & n_pwc == 1,
    DS_only      = has_DS     & n_pwc == 1,
    SS_only      = has_SS     & n_pwc == 1,
    VCBLDU_only  = has_VCBLDU & n_pwc == 1,
    HZ_only      = has_HZ     & n_pwc == 1,
    
    # Any & Multiple within the same report
    ANY_PWC_report      = n_pwc >= 1,
    MULTIPLE_PWC_report = n_pwc >= 2
  ) %>%
  group_by(station, date) %>%
  summarise(
    # counts of reports meeting each condition
    DU_only      = sum(DU_only,      na.rm = TRUE),
    DS_only      = sum(DS_only,      na.rm = TRUE),
    SS_only      = sum(SS_only,      na.rm = TRUE),
    VCBLDU_only  = sum(VCBLDU_only,  na.rm = TRUE),
    HZ_only      = sum(HZ_only,      na.rm = TRUE),
    
    ANY_PWC      = sum(ANY_PWC_report,      na.rm = TRUE),
    MULTIPLE_PWC = sum(MULTIPLE_PWC_report, na.rm = TRUE),
    
    # 0/1 daily flags (did it occur at least once today?)
    DU_only_flag      = as.integer(any(DU_only)),
    DS_only_flag      = as.integer(any(DS_only)),
    SS_only_flag      = as.integer(any(SS_only)),
    VCBLDU_only_flag  = as.integer(any(VCBLDU_only)),
    HZ_only_flag      = as.integer(any(HZ_only)),
    ANY_PWC_flag      = as.integer(any(ANY_PWC_report)),
    MULTIPLE_PWC_flag = as.integer(any(MULTIPLE_PWC_report)),
    .groups = "drop"
  ) %>%
  arrange(station, date)



#--------------------------------------   ------------------------------------------
# --- daily flags (add SS_day) ---
daily_flags5 <- dust_events_cv5 %>%
  mutate(
    date   = as.Date(local_time),
    txt_wx = paste(coalesce(as.character(wxcodes), ""),
                   coalesce(as.character(metar),   "")),
    
    has_DU     = str_detect(txt_wx, regex("\\bDU\\b",     TRUE)),
    has_BLDU   = str_detect(txt_wx, regex("\\bBLDU\\b",   TRUE)),
    has_VCBLDU = str_detect(txt_wx, regex("\\bVCBLDU\\b", TRUE)),
    has_DS     = str_detect(txt_wx, regex("\\bDS\\b",     TRUE)),
    has_SS     = str_detect(txt_wx, regex("\\bSS\\b",     TRUE)),
    has_HZ     = str_detect(txt_wx, regex("\\bHZ\\b",     TRUE))
  ) %>%
  mutate(
    has_BLDU_plain = has_BLDU & !has_VCBLDU,
    n_pwc = rowSums(cbind(has_DU, has_BLDU_plain, has_DS, has_SS, has_VCBLDU, has_HZ), na.rm = TRUE),
    ANY_PWC_report      = n_pwc >= 1,
    MULTIPLE_PWC_report = n_pwc >= 2
  ) %>%
  group_by(station, date) %>%
  summarise(
    DU_day        = as.integer(any(has_DU)),
    BLDU_day      = as.integer(any(has_BLDU_plain)),
    VCBLDU_day    = as.integer(any(has_VCBLDU)),
    DS_day        = as.integer(any(has_DS)),
    SS_day        = as.integer(any(has_SS)),      # << add this
    HZ_day        = as.integer(any(has_HZ)),
    ANY_day       = as.integer(any(ANY_PWC_report)),          # includes SS
    MULTIPLE_day  = as.integer(any(MULTIPLE_PWC_report)),
    .groups = "drop"
  )

# --- per-station totals (add SS_days) ---
station_totals_days5 <- daily_flags5 %>%
  group_by(station) %>%
  summarise(
    DU_days            = sum(DU_day),
    BLDU_days          = sum(BLDU_day),
    VCBLDU_days        = sum(VCBLDU_day),
    DS_days            = sum(DS_day),
    SS_days            = sum(SS_day),            # << add this
    HZ_days            = sum(HZ_day),
    ANY_days           = sum(ANY_day),
    MULTIPLE_PWC_days  = sum(MULTIPLE_day),
    .groups = "drop"
  ) %>%
  arrange(desc(ANY_days))


#---------------------------------------------------------start--------------------------------
library(dplyr)
library(stringr)
library(lubridate)

# --- Tunable thresholds for the midnight-continuation heuristic ---
POST_MIDNIGHT_CUTOFF_HOUR <- 1   # consider 00:00–01:59 as "after midnight"
CROSS_GAP_HOURS           <- 2   # max hours allowed across midnight to tag continuation

# ============ 1) Per-report flags ============
report_flags <- dust_events_cv5 %>%
  mutate(
    date   = as.Date(local_time),
    hour   = lubridate::hour(local_time),
    txt_wx = paste(coalesce(as.character(wxcodes), ""),
                   coalesce(as.character(metar),   "")),
    
    # Base detections (strict tokens, case-insensitive)
    has_DU     = str_detect(txt_wx, regex("\\bDU\\b",     TRUE)),
    has_BLDU   = str_detect(txt_wx, regex("\\bBLDU\\b",   TRUE)),
    has_VCBLDU = str_detect(txt_wx, regex("\\bVCBLDU\\b", TRUE)),
    has_DS     = str_detect(txt_wx, regex("\\bDS\\b",     TRUE)),
    has_SS     = str_detect(txt_wx, regex("\\bSS\\b",     TRUE)),
    has_HZ     = str_detect(txt_wx, regex("\\bHZ\\b",     TRUE))
  ) %>%
  # Treat VCBLDU separately so it doesn't also count as BLDU
  mutate(
    has_BLDU_plain = has_BLDU & !has_VCBLDU,
    
    # Count how many PWC codes appear in the SAME report
    n_pwc = rowSums(cbind(has_DU, has_BLDU_plain, has_DS, has_SS, has_VCBLDU, has_HZ), na.rm = TRUE),
    
    # "ONLY" = exactly one of the codes present in this report
    DU_only_report      = has_DU     & n_pwc == 1,
    DS_only_report      = has_DS     & n_pwc == 1,
    SS_only_report      = has_SS     & n_pwc == 1,
    BLDU_only_report    = has_BLDU_plain & n_pwc == 1,
    VCBLDU_only_report  = has_VCBLDU & n_pwc == 1,
    HZ_only_report      = has_HZ     & n_pwc == 1,
    
    ANY_PWC_report      = n_pwc >= 1,
    MULTIPLE_PWC_report = n_pwc >= 2
  )

# ============ 2) Continuation across midnight (per report, then per day) ============
pwc_seq <- report_flags %>%
  filter(ANY_PWC_report) %>%
  arrange(station, local_time) %>%
  group_by(station) %>%
  mutate(
    prev_time = lag(local_time),
    prev_date = as.Date(prev_time),
    time_diff_hours = as.numeric(difftime(local_time, prev_time, units = "hours")),
    # continued from previous day across midnight?
    entered_from_prev_report =
      !is.na(prev_time) &
      (as.Date(local_time) == prev_date + 1) &
      (hour(local_time) <= POST_MIDNIGHT_CUTOFF_HOUR) &
      (time_diff_hours <= CROSS_GAP_HOURS)
  ) %>%
  ungroup()

continuation_by_day <- pwc_seq %>%
  group_by(station, date = as.Date(local_time)) %>%
  summarise(entered_from_prev = as.integer(any(entered_from_prev_report, na.rm = TRUE)),
            .groups = "drop")

# ============ 3) Per-day table with ONLY/MULTIPLE counts & period tag ============
day_table <- report_flags %>%
  group_by(station, date) %>%
  summarise(
    DU_only       = sum(DU_only_report,      na.rm = TRUE),
    DS_only       = sum(DS_only_report,      na.rm = TRUE),
    SS_only       = sum(SS_only_report,      na.rm = TRUE),
    BLDU_only     = sum(BLDU_only_report,    na.rm = TRUE),
    VCBLDU_only   = sum(VCBLDU_only_report,  na.rm = TRUE),
    HZ_only       = sum(HZ_only_report,      na.rm = TRUE),
    MULTIPLE_PWC  = sum(MULTIPLE_PWC_report, na.rm = TRUE),
    
    # Period classification based on ANY PWC times that day
    morning_any   = any(ANY_PWC_report & hour < 12,  na.rm = TRUE),
    evening_any   = any(ANY_PWC_report & hour >= 12, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  ungroup() %>%
  mutate(
    period = case_when(
      morning_any & !evening_any ~ "Morning only",
      !morning_any & evening_any ~ "Evening only",
      morning_any & evening_any  ~ "Both",
      TRUE                       ~ "None"
    )
  ) %>%
  select(-morning_any, -evening_any) %>%
  left_join(continuation_by_day, by = c("station","date")) %>%
  mutate(entered_from_prev = if_else(is.na(entered_from_prev), 0L, entered_from_prev)) %>%
  arrange(station, date)

# >>> This is your detailed table per station/day <<<
#day_table
# Columns:
# station | date | DU_only | DS_only | SS_only | BLDU_only | VCBLDU_only | HZ_only | MULTIPLE_PWC | period | entered_from_prev

# ============ 4) Groupings you asked for ============

# (A) Group days by whether they entered from previous day
by_continuation <- day_table %>%
  mutate(
    DU_only_day      = as.integer(DU_only      > 0),
    DS_only_day      = as.integer(DS_only      > 0),
    SS_only_day      = as.integer(SS_only      > 0),
    BLDU_only_day    = as.integer(BLDU_only    > 0),
    VCBLDU_only_day  = as.integer(VCBLDU_only  > 0),
    HZ_only_day      = as.integer(HZ_only      > 0),
    MULTIPLE_PWC_day = as.integer(MULTIPLE_PWC > 0)
  ) %>%
  group_by(station, entered_from_prev) %>%
  summarise(
    DU_only_days      = sum(DU_only_day),
    DS_only_days      = sum(DS_only_day),
    SS_only_days      = sum(SS_only_day),
    BLDU_only_days    = sum(BLDU_only_day),
    VCBLDU_only_days  = sum(VCBLDU_only_day),
    HZ_only_days      = sum(HZ_only_day),
    MULTIPLE_PWC_days = sum(MULTIPLE_PWC_day),
    .groups = "drop"
  )



# (B) Group days by period (Morning/Evening/Both)
by_period <- day_table %>%
  mutate(
    DU_only_day      = as.integer(DU_only      > 0),
    DS_only_day      = as.integer(DS_only      > 0),
    SS_only_day      = as.integer(SS_only      > 0),
    BLDU_only_day    = as.integer(BLDU_only    > 0),
    VCBLDU_only_day  = as.integer(VCBLDU_only  > 0),
    HZ_only_day      = as.integer(HZ_only      > 0),
    MULTIPLE_PWC_day = as.integer(MULTIPLE_PWC > 0)
  ) %>%
  group_by(station, period) %>%
  summarise(
    DU_only_days      = sum(DU_only_day),
    DS_only_days      = sum(DS_only_day),
    SS_only_days      = sum(SS_only_day),
    BLDU_only_days    = sum(BLDU_only_day),
    VCBLDU_only_days  = sum(VCBLDU_only_day),
    HZ_only_days      = sum(HZ_only_day),
    MULTIPLE_PWC_days = sum(MULTIPLE_PWC_day),
    .groups = "drop"
  ) %>%
  arrange(station, period)



# --- Helper: filter one station’s day table (e.g., "BFL") ---
PTV_days <- day_table %>% filter(station == "PTV")

#--------------end-----------------------------------------------end-----------------------





#------------------------------------------------------
#------------------------------------------------------
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(purrr)

# ---------------- Config ----------------
MAX_GAP_HOURS <- 2  # split events when gap between PWC reports exceeds this

# ---------------- 0) Per-report flags + county map ----------------
base <- dust_events_cv5 %>%
  mutate(
    date   = as.Date(local_time),
    txt_wx = paste(coalesce(as.character(wxcodes), ""),
                   coalesce(as.character(metar),   "")),
    has_DU     = str_detect(txt_wx, regex("\\bDU\\b",     ignore_case = TRUE)),
    has_BLDU   = str_detect(txt_wx, regex("\\bBLDU\\b",   ignore_case = TRUE)),
    has_VCBLDU = str_detect(txt_wx, regex("\\bVCBLDU\\b", ignore_case = TRUE)),
    has_DS     = str_detect(txt_wx, regex("\\bDS\\b",     ignore_case = TRUE)),
    has_SS     = str_detect(txt_wx, regex("\\bSS\\b",     ignore_case = TRUE)),
    has_HZ     = str_detect(txt_wx, regex("\\bHZ\\b",     ignore_case = TRUE))
  ) %>%
  mutate(
    has_BLDU_plain = has_BLDU & !has_VCBLDU,
    n_pwc          = rowSums(cbind(has_DU, has_BLDU_plain, has_DS, has_SS, has_VCBLDU, has_HZ), na.rm = TRUE),
    ANY_PWC_report = n_pwc >= 1,
    
    county = case_when(
      station %in% c("OVE")        ~ "Butte",
      station %in% c("FAT")        ~ "Fresno",
      station ==  "BFL"            ~ "Kern",
      station ==  "HJO"            ~ "Kings",
      station ==  "MAE"            ~ "Madera",
      station ==  "MCE"            ~ "Merced",
      station %in% c("SAC","SMF")  ~ "Sacramento",
      station ==  "RDD"            ~ "Shasta",
      station ==  "BAB"            ~ "Sutter",
      station ==  "RBL"            ~ "Tehama",
      station %in% c("PTV","VIS")  ~ "Tulare",
      station ==  "MYV"            ~ "Yuba",
      station ==  "SCK"            ~ "San Joaquin",
      TRUE                         ~ NA_character_
    )
  )

# ---------------- Presence by day (station) ----------------
dust_presence <- base %>%
  group_by(station, date) %>%
  summarise(dust_day = as.integer(any(ANY_PWC_report, na.rm = TRUE)), .groups = "drop")

# ---------------- 1) Segment events per station ----------------
events <- base %>%
  filter(ANY_PWC_report) %>%
  arrange(station, local_time) %>%
  group_by(station) %>%
  mutate(
    dt_hours  = as.numeric(difftime(local_time, lag(local_time), units = "hours")),
    new_event = is.na(dt_hours) | dt_hours > MAX_GAP_HOURS,
    event_id  = cumsum(new_event)
  ) %>%
  ungroup() %>%
  group_by(station, event_id) %>%
  summarise(
    start_time = min(local_time),
    end_time   = max(local_time),
    .groups = "drop"
  ) %>%
  mutate(
    start_date       = as.Date(start_time),
    end_date         = as.Date(end_time),
    crosses_midnight = start_date != end_date,
    start_period     = if_else(hour(start_time) < 12, "Morning", "Evening")
  )

# ---------------- 2) Expand station events across dates ----------------
event_cover <- events %>%
  mutate(date_seq = map2(start_date, end_date, ~ seq(.x, .y, by = "day"))) %>%
  unnest(date_seq) %>%
  rename(date = date_seq) %>%
  mutate(
    started_today   = (date == start_date),
    enters_next_day = started_today & crosses_midnight
  )

# ---------------- 3) FIRST TABLE: per station–date ----------------
day_events_table <- dust_presence %>%
  left_join(
    event_cover %>%
      group_by(station, date) %>%
      summarise(
        morning_event       = as.integer(any(started_today & start_period == "Morning")),
        evening_event       = as.integer(any(started_today & start_period == "Evening")),
        continued_from_prev = as.integer(any(!started_today)),
        enters_next_day     = as.integer(any(enters_next_day)),
        .groups = "drop"
      ),
    by = c("station","date")
  ) %>%
  mutate(
    morning_event       = replace_na(morning_event, 0L),
    evening_event       = replace_na(evening_event, 0L),
    continued_from_prev = replace_na(continued_from_prev, 0L),
    enters_next_day     = replace_na(enters_next_day, 0L),
    events_today        = morning_event + evening_event
  ) %>%
  arrange(station, date)

# ---------------- 4) Station totals ----------------
station_events_total <- day_events_table %>%
  group_by(station) %>%
  summarise(
    total_dust_days       = sum(dust_day,       na.rm = TRUE),
    total_events          = sum(events_today,   na.rm = TRUE),
    morning_events_total  = sum(morning_event,  na.rm = TRUE),
    evening_events_total  = sum(evening_event,  na.rm = TRUE),
    cross_midnight_events = sum(enters_next_day,na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_events))

# =================== COUNTY DE-DUP BY DAY & PERIOD ===================
# Count one event per county per date per period (Morning/Evening).
station_to_county <- base %>%
  distinct(station, county) %>%
  filter(!is.na(county))

day_events_with_county <- day_events_table %>%
  left_join(station_to_county, by = "station") %>%
  filter(!is.na(county))

county_day_events_table <- day_events_with_county %>%
  group_by(county, date) %>%
  summarise(
    dust_day_county     = as.integer(any(dust_day == 1)),
    morning_event       = as.integer(any(morning_event == 1)),  # ≤1 per county-day
    evening_event       = as.integer(any(evening_event == 1)),  # ≤1 per county-day
    events_today        = morning_event + evening_event,        # 0..2
    continued_from_prev = as.integer(any(continued_from_prev == 1)),
    enters_next_day     = as.integer(any(enters_next_day == 1)),
    .groups = "drop"
  ) %>%
  arrange(county, date)

county_events_total <- county_day_events_table %>%
  group_by(county) %>%
  summarise(
    total_dust_days      = sum(dust_day_county, na.rm = TRUE),
    total_events         = sum(events_today,     na.rm = TRUE),
    morning_events_total = sum(morning_event,    na.rm = TRUE),
    evening_events_total = sum(evening_event,    na.rm = TRUE),
    cross_midnight_events= sum(enters_next_day,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_events))

# (Optional) grand totals across all counties
all_counties_totals <- county_day_events_table %>%
  summarise(
    total_dust_days      = sum(dust_day_county, na.rm = TRUE),
    total_events         = sum(events_today,     na.rm = TRUE),
    morning_events_total = sum(morning_event,    na.rm = TRUE),
    evening_events_total = sum(evening_event,    na.rm = TRUE),
    cross_midnight_events= sum(enters_next_day,  na.rm = TRUE)
  )

# ----- Outputs -----
day_events_table
station_events_total
county_day_events_table
county_events_total
all_counties_totals



#------------------------------------days in multiple county------ptv vis smf sac


# --- station -> county map from your 'base' ---
station_to_county <- base %>%
  dplyr::distinct(station, county) %>%
  dplyr::filter(!is.na(county))

# --- pick the two counties where de-dup applies ---
target_counties <- c("Sacramento","Tulare")  # (SAC+SMF) and (PTV+VIS)

# --- for each station/day/period, keep the earliest start time that day ---
station_period_starts <- events %>%
  dplyr::transmute(
    station,
    start_time,
    date   = as.Date(start_time),
    period = dplyr::if_else(lubridate::hour(start_time) < 12, "Morning", "Evening")
  ) %>%
  dplyr::group_by(station, date, period) %>%
  dplyr::summarise(start_time = min(start_time), .groups = "drop") %>%
  dplyr::left_join(station_to_county, by = "station") %>%
  dplyr::filter(county %in% target_counties)

# --- (A) Compact table: one row per county/date/period where ≥2 stations started an event ---
county_period_collisions <- station_period_starts %>%
  dplyr::group_by(county, date, period) %>%
  dplyr::filter(dplyr::n() >= 2) %>%
  dplyr::summarise(
    stations     = paste(sort(unique(station)), collapse = ", "),
    start_times  = paste(format(start_time, "%H:%M"), collapse = ", "),
    first_start  = min(start_time),
    last_start   = max(start_time),
    .groups = "drop"
  ) %>%
  dplyr::arrange(county, date, period)

# --- (B) Wide view: columns SAC/SMF (and PTV/VIS) with their times ---
county_period_wide <- station_period_starts %>%
  dplyr::semi_join(county_period_collisions, by = c("county","date","period")) %>%
  dplyr::mutate(start_time_local = format(start_time, "%Y-%m-%d %H:%M")) %>%
  dplyr::select(county, date, period, station, start_time_local) %>%
  tidyr::pivot_wider(names_from = station, values_from = start_time_local) %>%
  dplyr::arrange(county, date, period)

# --- Inspect tables ---
county_period_collisions
county_period_wide










#----------------------------------------

library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(forcats)

# ---------------- CONFIG ----------------
MAX_GAP_HOURS <- 2
DUR_LABELS <- c("≤1h","2h","3h","4h","5h","6h","7h","8h","9h","≥10h")
DUR_BREAKS <- c(-Inf, 1, 2, 3, 4, 5, 6, 7, 8, 9, Inf)

# ---------------- 1) Per-report PWC flags (station level) ----------------
base <- dust_events_cv5 %>%
  mutate(
    txt_wx = paste(coalesce(as.character(wxcodes), ""),
                   coalesce(as.character(metar),   "")),
    has_DU     = str_detect(txt_wx, regex("\\bDU\\b",     TRUE)),
    has_BLDU   = str_detect(txt_wx, regex("\\bBLDU\\b",   TRUE)),
    has_VCBLDU = str_detect(txt_wx, regex("\\bVCBLDU\\b", TRUE)),
    has_DS     = str_detect(txt_wx, regex("\\bDS\\b",     TRUE)),
    has_SS     = str_detect(txt_wx, regex("\\bSS\\b",     TRUE)),
    has_HZ     = str_detect(txt_wx, regex("\\bHZ\\b",     TRUE))
  ) %>%
  mutate(
    has_BLDU_plain = has_BLDU & !has_VCBLDU,
    n_pwc          = rowSums(cbind(has_DU, has_BLDU_plain, has_DS, has_SS, has_VCBLDU, has_HZ), na.rm = TRUE),
    ANY_PWC_report = n_pwc >= 1
  ) %>%
  filter(ANY_PWC_report) %>%
  arrange(station, local_time)

# ---------------- 2) Segment EVENTS per station ----------------
events <- base %>%
  group_by(station) %>%
  mutate(
    dt_hours  = as.numeric(difftime(local_time, lag(local_time), units = "hours")),
    new_event = is.na(dt_hours) | dt_hours > MAX_GAP_HOURS,
    event_id  = cumsum(new_event)
  ) %>%
  ungroup() %>%
  group_by(station, event_id) %>%
  summarise(
    start_time = min(local_time),
    end_time   = max(local_time),
    .groups = "drop"
  ) %>%
  mutate(
    start_date   = as.Date(start_time),
    end_date     = as.Date(end_time),
    start_year   = year(start_time),
    start_month  = factor(month.abb[month(start_time)], levels = month.abb, ordered = TRUE),
    start_hour   = hour(start_time),
    duration_hr  = as.numeric(difftime(end_time, start_time, units = "hours")),
    duration_bin = cut(duration_hr, breaks = DUR_BREAKS, labels = DUR_LABELS, right = TRUE, ordered_result = TRUE),
    start_period = if_else(start_hour < 12, "Morning", "Evening")
  )

# ---------------- 3) Build the per-station day table (max 2 events/day) ----------------
# presence by day (dust day)
dust_presence <- base %>%
  mutate(date = as.Date(local_time)) %>%
  group_by(station, date) %>%
  summarise(dust_day = as.integer(any(ANY_PWC_report)), .groups = "drop")

# cover each event across its days to get continued/enter-next info
event_cover <- events %>%
  mutate(date_seq = Map(seq, start_date, end_date, by = "day")) %>%
  tidyr::unnest(date_seq) %>%
  rename(date = date_seq) %>%
  mutate(
    started_today   = (date == start_date),
    enters_next_day = started_today & (start_date != end_date)
  )

day_events_table <- dust_presence %>%
  left_join(
    event_cover %>%
      group_by(station, date) %>%
      summarise(
        morning_event       = as.integer(any(started_today & start_period == "Morning")),
        evening_event       = as.integer(any(started_today & start_period == "Evening")),
        continued_from_prev = as.integer(any(!started_today)),
        enters_next_day     = as.integer(any(enters_next_day)),
        .groups = "drop"
      ),
    by = c("station","date")
  ) %>%
  mutate(
    across(c(morning_event, evening_event, continued_from_prev, enters_next_day),
           ~ replace_na(.x, 0L)),
    events_today = morning_event + evening_event
  ) %>%
  arrange(station, date)

# ---------------- 4) Station totals (unchanged) ----------------
station_events_total <- day_events_table %>%
  group_by(station) %>%
  summarise(
    total_dust_days       = sum(dust_day),
    total_events          = sum(events_today),
    morning_events_total  = sum(morning_event),
    evening_events_total  = sum(evening_event),
    cross_midnight_events = sum(enters_next_day),
    .groups = "drop"
  ) %>%
  arrange(desc(total_events))

# ---------------- 5) COUNTY DE-DUP BY DAY & PERIOD ----------------
# Map station -> county
station_to_county <- tibble::tibble(
  station = c("OVE","FAT","BFL","HJO","MAE","MCE","SAC","SMF","RDD","BAB","RBL","PTV","VIS","MYV","SCK"),
  county  = c("Butte","Fresno","Kern","Kings","Madera","Merced","Sacramento","Sacramento","Shasta",
              "Sutter","Tehama","Tulare","Tulare","Yuba","San Joaquin")
)

# Merge county into station-day table
day_events_with_county <- day_events_table %>%
  left_join(station_to_county, by = "station") %>%
  filter(!is.na(county))

# >>> County de-dup logic: for each county+date, OR across stations within Morning/Evening
county_day_events_table <- day_events_with_county %>%
  group_by(county, date) %>%
  summarise(
    dust_day_county   = as.integer(any(dust_day == 1)),              # any station had dust that day
    morning_event     = as.integer(any(morning_event == 1)),          # morning counted once per county-day
    evening_event     = as.integer(any(evening_event == 1)),          # evening counted once per county-day
    events_today      = morning_event + evening_event,                # 0..2
    continued_from_prev = as.integer(any(continued_from_prev == 1)),  # optional metadata
    enters_next_day     = as.integer(any(enters_next_day == 1)),      # optional metadata
    .groups = "drop"
  ) %>%
  arrange(county, date)

# County totals with de-dup (each county-day Morning/Evening counts once)
county_events_total <- county_day_events_table %>%
  group_by(county) %>%
  summarise(
    total_dust_days      = sum(dust_day_county),
    total_events         = sum(events_today),
    morning_events_total = sum(morning_event),
    evening_events_total = sum(evening_event),
    cross_midnight_events= sum(enters_next_day),
    .groups = "drop"
  ) %>%
  arrange(desc(total_events))

# ---------------- 6) (Optional) overall totals across all counties ----------------
all_counties_totals <- county_day_events_table %>%
  summarise(
    total_dust_days      = sum(dust_day_county),
    total_events         = sum(events_today),
    morning_events_total = sum(morning_event),
    evening_events_total = sum(evening_event),
    cross_midnight_events= sum(enters_next_day)
  )

# -------- Outputs --------
day_events_table                     # station-day
station_events_total                 # station totals
county_day_events_table              # county-day (dedup by Morning/Evening)
county_events_total                  # county totals
all_counties_totals                  # one-row grand totals



#----------------------------------------------------------


# ================== COUNTY-DE-DUPPED EVENT DISTRIBUTIONS ==================


# station -> county (already in `base`)
station_to_county <- base %>%
  dplyr::distinct(station, county) %>%
  dplyr::filter(!is.na(county))

# attach county to station events
events_with_county <- events %>%
  dplyr::left_join(station_to_county, by = "station") %>%
  dplyr::filter(!is.na(county))

# collapse to ONE county event per (county, start_date, start_period)
county_events_dedup <- events_with_county %>%
  dplyr::group_by(county, start_date, start_period) %>%
  dplyr::summarise(
    start_time = min(start_time, na.rm = TRUE),
    end_time   = max(end_time,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    start_year   = lubridate::year(start_time),
    start_month  = factor(month.abb[lubridate::month(start_time)], levels = month.abb, ordered = TRUE),
    start_hour   = lubridate::hour(start_time),
    duration_hr  = as.numeric(difftime(end_time, start_time, units = "hours")),
    duration_bin = cut(duration_hr,
                       breaks = DUR_BREAKS,
                       labels = DUR_LABELS,
                       right = TRUE,
                       ordered_result = TRUE)
  )

# -------- Totals across ALL counties (dedup already applied) --------
events_by_year_total <- county_events_dedup %>%
  dplyr::count(start_year, name = "events") %>%
  dplyr::arrange(start_year)

events_by_month_total <- county_events_dedup %>%
  dplyr::count(start_month, name = "events") %>%
  tidyr::complete(start_month = factor(month.abb, levels = month.abb, ordered = TRUE),
                  fill = list(events = 0)) %>%
  dplyr::arrange(start_month)

events_by_hour_total <- county_events_dedup %>%
  dplyr::count(start_hour, name = "events") %>%
  tidyr::complete(start_hour = 0:23, fill = list(events = 0)) %>%
  dplyr::arrange(start_hour)

events_by_duration_total <- county_events_dedup %>%
  dplyr::count(duration_bin, name = "events") %>%
  tidyr::complete(duration_bin = factor(DUR_LABELS, levels = DUR_LABELS, ordered = TRUE),
                  fill = list(events = 0)) %>%
  dplyr::arrange(duration_bin)

# -------- Bundle (exact names you used) --------
tables_total <- list(
  yearly   = events_by_year_total,
  monthly  = events_by_month_total,
  diurnal  = events_by_hour_total,
  duration = events_by_duration_total
)


#--------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(maps)
library(cowplot)
library(gridExtra)
library(grid) 
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
library(ggrepel)
library(ggspatial)


# 1) Read & reproject the boundary
cv_boundary <- st_read("central_valley_alluvial_boundary_usgs.shp", quiet = TRUE) %>%
  st_zm(drop=TRUE) %>%
  st_transform(4326)


year_totals5=tables_total$yearly
colnames(year_totals5)=c('year','count')

month_totals5=tables_total$monthly
colnames(month_totals5)=c('month','count')

dur_counts5=tables_total$duration
colnames(dur_counts5)=c('duration_group','count')

diurnal_totals5=tables_total$diurnal
colnames(diurnal_totals5)=c('hour','count')

#-------------------------------------------------------------------
library(scales)

# make sure the factor order is what you want
dur_counts5 <- dur_counts5 %>%
  mutate(duration_group = factor(duration_group,
                                 levels = c("≤1h","2h","3h","4h","5h","6h","7h","8h","9h","≥10h"),
                                 ordered = TRUE)) %>%
  arrange(duration_group) %>%
  mutate(
    pct     = count / sum(count),
    pct_lab = percent(pct, accuracy = 0.1)  # e.g., "74.5%"
  )



#––– Panel C: Duration + inset –––
ymax_c <- max(dur_counts5$count) * 1.15  # leave headroom for labels
#-----------------------------------------------------------------------


#––– 5) Prepare inset map –––
cv_counties <- map_data("county", region = "california") %>%
  filter(subregion %in% c("shasta","tehama","glenn","butte","colusa",
                          "sutter","yuba","yolo","sacramento","san joaquin",
                          "stanislaus","merced","madera","fresno","kings",
                          "tulare","kern"))
ca_state <- map_data("state", region = "california")

map_insert1=ggplot() +
  # a) Central Valley alluvial boundary base (your sf layer)
  geom_sf(data = cv_boundary,
          fill  = "#f0f0f0",
          color = "black",
          linewidth = 0.6) +
  # b) county subdivisions (your fortified polygon layer)
  geom_polygon(data = cv_counties,
               aes(x = long, y = lat, group = group),
               fill  = NA,
               color = "red",
               linewidth = 0.1)+
  theme_void() 

map_inset <- ggplot() +
  geom_polygon(data = ca_state, aes(long, lat, group = group),
               fill="grey90", color="black") +
  geom_polygon(data = cv_counties, aes(long, lat, group = group),
               fill="lightblue", color="blue") +
  coord_quickmap() +
  theme_void()

#––– Theme helper –––
clean_theme <- theme_classic(base_size = 24) +
  theme(
    panel.border   = element_blank(),
    axis.line      = element_blank(),
    axis.ticks     = element_blank()
  )

#––– Panel A: Yearly w/ trend –––
pA <- ggplot(year_totals5, aes(year, count)) +
  geom_col(fill="black") +
  geom_smooth(method="lm", se=FALSE, linetype="dashed", color="black") +
  ylim(0,120)+
  labs(title="(a) Yearly Dust Event (2005–2024)",
       x="Year", y=NULL) +
  clean_theme +
  theme(plot.title = element_text(hjust = 0))

#––– Panel B: Monthly –––
pB <- ggplot(month_totals5, aes(month, count)) +
  geom_col(fill="black") +ylim(0,220)+
  labs(title="(b) Monthly Dust Event",
       x="Month", y=NULL) +
  clean_theme +
  theme(plot.title = element_text(hjust = 0))

#––– Panel C: Duration + inset –––
pC_base <- ggplot(dur_counts5, aes(duration_group, count)) +
  geom_col(fill="black") +
  labs(title="(C) Event Count by Duration",
       x="Duration (Hours)", y=NULL) +
  clean_theme +
  theme(plot.title = element_text(hjust = 0))

pC <- ggdraw(pC_base) +
  draw_plot(map_insert1, x=0.55, y=0.5, width=0.4, height=0.4)

#--------------------------------------------------------------------
pC_base <- ggplot(dur_counts5, aes(duration_group, count)) +
  geom_col(fill="black") +
  geom_text(aes(label = pct_lab), vjust = -0.35, size = 5) +
  coord_cartesian(ylim = c(0, ymax_c)) +
  labs(title="(c) Event Count by Duration",
       x="Duration (Hours)", y=NULL) +
  clean_theme +
  theme(plot.title = element_text(hjust = 0))

pC <- ggdraw(pC_base) +
  draw_plot(map_insert1, x=0.55, y=0.5, width=0.4, height=0.4)
#---------------------------------------------------------------------

#––– Panel D: Diurnal –––
pD <- ggplot(diurnal_totals5, aes(hour, count)) +
  geom_col(fill="black", width=1) +
  scale_x_continuous(breaks=seq(0,24,2),
                     labels=sprintf("%02d", seq(0,24,2))) +
  coord_cartesian(xlim=c(0,23)) +ylim(0,150)+
  labs(title="(d) Diurnal Start-Time Distribution",
       x="Local Time (Hours)", y=NULL) +
  clean_theme +
  theme(plot.title = element_text(hjust = 0))

#––– 6) Arrange and add single Y label –––
png("devvvvall_plots9.png", width=1850, height=950)
grid.arrange(
  pA, pB,
  pC, pD,
  ncol = 2, nrow = 2,
  left = textGrob("Total Dust Events", rot = 90,
                  gp = gpar(fontsize = 26, fontface="bold"))
)
dev.off()


#-----------------plot duration fig s2--------------------------------------------------------
# ------------------------------------------------------
# libs
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(forcats)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(sf)

# ------------------------------------------------------
# CONFIG
MAX_GAP_HOURS <- 2
DUR_LABELS <- c(as.character(1:9), ">10")
DUR_BREAKS <- c(-Inf, 1, 2, 3, 4, 5, 6, 7, 8, 9, Inf)  # ≤1, 2..9, ≥10

# ------------------------------------------------------
# 1) Build per-report PWC flags from dust_events_cv5 (no wind filter here)
base_cv <- dust_events_cv5 %>%
  mutate(
    txt_wx = paste(coalesce(as.character(wxcodes), ""),
                   coalesce(as.character(metar), "")),
    has_DU     = str_detect(txt_wx, regex("\\bDU\\b",     TRUE)),
    has_BLDU   = str_detect(txt_wx, regex("\\bBLDU\\b",   TRUE)),
    has_VCBLDU = str_detect(txt_wx, regex("\\bVCBLDU\\b", TRUE)),
    has_DS     = str_detect(txt_wx, regex("\\bDS\\b",     TRUE)),
    has_SS     = str_detect(txt_wx, regex("\\bSS\\b",     TRUE)),
    has_HZ     = str_detect(txt_wx, regex("\\bHZ\\b",     TRUE))
  ) %>%
  mutate(
    has_BLDU_plain = has_BLDU & !has_VCBLDU,
    n_pwc          = rowSums(cbind(has_DU, has_BLDU_plain, has_DS, has_SS, has_VCBLDU, has_HZ), na.rm = TRUE),
    ANY_PWC_report = n_pwc >= 1
  ) %>%
  filter(ANY_PWC_report) %>%
  arrange(station, local_time)

# ------------------------------------------------------
# 2) Segment events per station, compute durations & bins
dust_event_duration <- base_cv %>%
  group_by(station) %>%
  mutate(
    dt_hours  = as.numeric(difftime(local_time, lag(local_time), units = "hours")),
    new_event = is.na(dt_hours) | dt_hours > MAX_GAP_HOURS,
    event_id  = cumsum(new_event)
  ) %>%
  ungroup() %>%
  group_by(station, event_id) %>%
  summarise(
    start_time = min(local_time),
    end_time   = max(local_time),
    .groups    = "drop"
  ) %>%
  mutate(
    start_date   = as.Date(start_time),
    start_year   = year(start_time),
    start_month  = factor(month.abb[month(start_time)], levels = month.abb, ordered = TRUE),
    start_hour   = hour(start_time),
    start_period = if_else(start_hour < 12, "Morning", "Evening"),
    duration_hr  = as.numeric(difftime(end_time, start_time, units = "hours")),
    duration_group = cut(duration_hr, breaks = DUR_BREAKS, labels = DUR_LABELS,
                         right = TRUE, ordered_result = TRUE)
  ) %>%
  arrange(station, start_time)

# ------------------------------------------------------
# --- Station coordinates (order north → south) ---
coords_txt <- "
station,lat,lon
BFL,35.43440,-119.0542
BAB,39.13609,-121.4366
FAT,36.78000,-119.7194
HJO,36.31139,-119.6232
MAE,36.98486,-120.1107
MYV,39.10203,-121.5688
MCE,37.28603,-120.5179
OVE,39.49000,-121.6200
PTV,36.02732,-119.0629
RBL,40.15190,-122.2536
RDD,40.50900,-122.2934
SAC,38.50690,-121.4950
SMF,38.69542,-121.5908
SCK,37.89417,-121.2383
VIS,36.31867,-119.3929
"
df_coords <- read.csv(text = coords_txt, stringsAsFactors = FALSE) %>%
  arrange(desc(lat))                     # north → south
station_ids <- intersect(df_coords$station, unique(dust_event_duration$station))
df_coords   <- df_coords %>% filter(station %in% station_ids)
station_ids <- df_coords$station

# --- Central Valley boundary (for the map inset) ---
library(sf); library(ggplot2); library(cowplot); library(gridExtra); library(forcats); library(tidyr); library(dplyr)

cv_boundary <- st_read("central_valley_alluvial_boundary_usgs.shp", quiet = TRUE) |>
  st_zm(drop = TRUE) |>
  st_transform(4326)

bb   <- sf::st_bbox(cv_boundary)
pad  <- 0.25
xlim <- c(bb["xmin"] - pad, bb["xmax"] + pad)
ylim <- c(bb["ymin"] - pad, bb["ymax"] + pad)

# one inset map per station
inset_maps <- setNames(vector("list", length(station_ids)), station_ids)
for (i in seq_len(nrow(df_coords))) {
  sc <- df_coords[i, ]
  inset_maps[[sc$station]] <-
    ggplot() +
    geom_sf(data = cv_boundary, fill = "#f0f0f0", color = "black", linewidth = 0.6) +
    geom_point(data = sc, aes(lon, lat), color = "red", size = 3, inherit.aes = FALSE) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void()
}

# --- Build two-bin (≤1 hr, >1 hr) % bars per station ---
# duration_group == "1" corresponds to (0, 1] because we used right-closed bins earlier.
plots_2bin <- list()

for (stn in station_ids) {
  stn_summary <- dust_event_duration %>%
    filter(station == stn) %>%
    mutate(
      bin = if_else(as.character(duration_group) == "1", "≤ 1 hr", "> 1 hr"),
      bin = factor(bin, levels = c("≤ 1 hr", "> 1 hr"), ordered = TRUE)
    ) %>%
    count(bin, name = "count") %>%
    complete(bin, fill = list(count = 0)) %>%
    mutate(
      total  = sum(count, na.rm = TRUE),
      pct    = ifelse(total > 0, 100 * count / total, 0),
      pct_lab = sprintf("%.1f%%", pct)
    )
  
  p_main <- ggplot(stn_summary, aes(x = bin, y = pct)) +
    geom_col(fill = "black", width = 0.8) +
    geom_text(aes(label = pct_lab), vjust = -0.35, size = 16) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = stn, x = NULL, y = "Count(%)") +
    theme_minimal(base_size = 26) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 26, face = "bold"),
      axis.text.x = element_text(size = 26),
      axis.text.y = element_text(size = 26)
    )
  
  inset <- inset_maps[[stn]]
  if (is.null(inset)) inset <- ggplot() + theme_void()
  
  plots_2bin[[stn]] <- cowplot::ggdraw() +
    cowplot::draw_plot(p_main) +
    cowplot::draw_plot(inset, x = 0.68, y = 0.60, width = 0.4, height = 0.4)
}

# --- Arrange grid (north → south) and save ---
n <- length(plots_2bin)
ncol <- 4
nrow <- ceiling(n / ncol)

png("station_duration_two_bins.png", width = 3000, height = 2000)
gridExtra::grid.arrange(grobs = plots_2bin, nrow = nrow, ncol = ncol)
dev.off()



#--------------------------------------------------------------









# ── PACKAGES ──────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(viridis)
  library(ggrepel)
  library(maps)
  library(patchwork)
  library(lubridate)
  library(grid)   # unit()
})

# ── INPUTS ────────────────────────────────────────────────────────
# Data frame `dust_events_cv5` with: station, local_time (POSIXct), wxcodes, metar
CV_SHP <- "central_valley_alluvial_boundary_usgs.shp"

# ── STATION COORDS (your list) ───────────────────────────────────
stations_df <- tibble::tribble(
  ~station, ~lat,      ~lon,
  "BFL",    35.43440, -119.0542,
  "BAB",    39.13609, -121.4366,
  "FAT",    36.78000, -119.7194,
  "HJO",    36.31139, -119.6232,
  "MAE",    36.98486, -120.1107,
  "MYV",    39.10203, -121.5688,
  "MCE",    37.28603, -120.5179,
  "OVE",    39.49000, -121.6200,
  "RBL",    40.15190, -122.2536,
  "RDD",    40.50900, -122.2934,
  "SAC",    38.50690, -121.4950,
  "SMF",    38.69542, -121.5908,
  "SCK",    37.89417, -121.2383,
  "VIS",    36.31867, -119.3929,
  "PTV",    36.02732, -119.0629
)

# ── 1) Parse codes (normalize VC BLDU -> VCBLDU) & collapse to DAILY PRESENCE ──
PWC_CODES <- c("BLDU","DU","DS","SS","HZ","VCBLDU")

norm_txt <- function(x) {
  x %>% stringr::str_to_upper() %>% stringr::str_replace_all("VC\\s*BLDU", "VCBLDU")
}
rx <- regex(paste0("\\b(", paste(PWC_CODES, collapse="|"), ")\\b"), ignore_case = TRUE)

per_report <- dust_events_cv5 %>%
  mutate(
    date    = as.Date(local_time),
    wxcodes = coalesce(as.character(wxcodes), ""),
    metar   = coalesce(as.character(metar),   ""),
    alltxt  = norm_txt(paste(wxcodes, metar, sep = " "))
  ) %>%
  mutate(
    codes_list  = str_extract_all(alltxt, rx)              %>% map(~ unique(toupper(.x))),
    n_codes_row = map_int(codes_list, length),
    multi_row   = n_codes_row >= 2
  ) %>%
  select(station, date, codes_list, multi_row)

# Daily presence booleans
daily_presence <- per_report %>%
  unnest(codes_list, keep_empty = TRUE) %>%
  filter(!is.na(codes_list)) %>%
  mutate(code = codes_list) %>%
  distinct(station, date, code) %>%
  pivot_wider(names_from = code, values_from = code, values_fn = length, values_fill = 0)

# Ensure all code columns exist
for (cc in PWC_CODES) if (!cc %in% names(daily_presence)) daily_presence[[cc]] <- 0L

daily_presence <- daily_presence %>%
  transmute(
    station, date,
    DU_day     = DU      > 0,
    BLDU_day   = BLDU    > 0,
    VCBLDU_day = VCBLDU  > 0,
    DS_day     = DS      > 0,
    SS_day     = SS      > 0,
    HZ_day     = HZ      > 0
  ) %>%
  # treat BLDU as "plain BLDU" (exclude VCBLDU)
  mutate(BLDU_day = BLDU_day & !VCBLDU_day)

# Multiple-PWC day (choose definition)
# Option A (default here): ≥2 DISTINCT codes ANYWHERE in the same day
multi_day_anywhere <- per_report %>%
  mutate(code_set = map(codes_list, ~ setdiff(.x, character(0)))) %>%
  group_by(station, date) %>%
  summarise(n_distinct_codes = length(unique(unlist(code_set))), .groups = "drop") %>%
  mutate(MULTI_day = n_distinct_codes >= 2)

# Option B: at least one report that had ≥2 codes simultaneously
multi_day_same_report <- per_report %>%
  group_by(station, date) %>%
  summarise(MULTI_same_report = any(multi_row), .groups = "drop")

# Merge presence + multi-day flags (use Option A by default)
daily_flags <- stations_df %>% select(station) %>%
  left_join(daily_presence,        by = "station") %>%
  left_join(multi_day_anywhere,    by = c("station","date")) %>%
  left_join(multi_day_same_report, by = c("station","date")) %>%
  mutate(
    across(c(DU_day,BLDU_day,VCBLDU_day,DS_day,SS_day,HZ_day), ~ replace_na(.x, FALSE)),
    MULTI_day         = replace_na(MULTI_day, FALSE),          # across-day (default)
    MULTI_same_report = replace_na(MULTI_same_report, FALSE)
  ) %>%
  mutate(ANY_dust_day = DU_day | BLDU_day | VCBLDU_day | DS_day | SS_day | HZ_day)

# ── 2) Per-station COUNTS (presence) ─────────────────────────────
counts_presence <- daily_flags %>%
  group_by(station) %>%
  summarise(
    total_days     = sum(ANY_dust_day, na.rm = TRUE),
    DU_days        = sum(DU_day,       na.rm = TRUE),
    BLDU_days      = sum(BLDU_day,     na.rm = TRUE),     # plain BLDU (excludes VCBLDU)
    DS_days        = sum(DS_day,       na.rm = TRUE),
    SS_days        = sum(SS_day,       na.rm = TRUE),
    HZ_days        = sum(HZ_day,       na.rm = TRUE),
    MULTI_PWC_days = sum(MULTI_day,    na.rm = TRUE),     # across-day definition
    # If you want the same-report metric instead, use:
    # MULTI_PWC_days = sum(MULTI_same_report, na.rm = TRUE),
    .groups = "drop"
  )

# ── 3) Geometry (CV boundary + counties) ─────────────────────────
cv_boundary <- st_read(CV_SHP, quiet = TRUE) %>%
  st_zm(drop = TRUE) %>%
  st_transform(4326)

central_valley_counties <- c(
  "shasta","tehama","glenn","butte","colusa",
  "sutter","yuba","yolo","sacramento","san joaquin",
  "stanislaus","merced","madera","fresno","kings",
  "tulare","kern"
)
ca_counties <- map_data("county", "california") %>%
  filter(subregion %in% central_valley_counties)

# ── 4) Map helper with safe color scale when vmax==0 ─────────────
make_map_from_vec <- function(vals_df, value_col, title_key,
                              option = "plasma", trans = "sqrt",
                              barheight_pt = 160) {
  
  dat <- stations_df %>%
    left_join(vals_df, by = "station") %>%
    mutate(value = !!rlang::sym(value_col), value = replace_na(value, 0L))
  
  vmax <- max(dat$value, na.rm = TRUE)
  vmax_adj <- if (is.finite(vmax) && vmax > 0) vmax else 1  # avoid 0/Inf warnings
  
  pos_df  <- filter(dat, value > 0)
  zero_df <- filter(dat, value == 0)
  sf_pos  <- st_as_sf(pos_df,  coords = c("lon","lat"), crs = 4326)
  sf_zero <- st_as_sf(zero_df, coords = c("lon","lat"), crs = 4326)
  
  p <- ggplot() +
    geom_sf(data = cv_boundary, fill = "#f0f0f0", color = "black", linewidth = 0.5) +
    geom_polygon(data = ca_counties, aes(long, lat, group = group),
                 fill = NA, color = "white", linewidth = 0.3) +
    geom_sf(data = sf_zero, shape = 21, fill = "white", color = "grey30", size = 6, linewidth = 0.6) +
    geom_sf(data = sf_pos,  aes(fill = value), shape = 21, color = "black", alpha = 0.9, size = 6) +
    scale_fill_viridis_c(
      name   = "Days",
      option = option, begin = 0.15, end = 0.95,
      limits = c(0, vmax_adj), oob = scales::squish,
      trans  = trans
    ) +
    geom_text_repel(data = dat, aes(lon, lat, label = station),
                    size = 3, segment.size = 0.2) +
    coord_sf(
      xlim = st_bbox(cv_boundary)[c("xmin","xmax")],
      ylim = st_bbox(cv_boundary)[c("ymin","ymax")],
      expand = FALSE
    ) +
    labs(title = title_key) +
    theme_void(base_size = 12) +
    guides(fill = guide_colorbar(title.position = "top",
                                 direction = "vertical",
                                 barheight = grid::unit(barheight_pt, "pt"),
                                 barwidth  = grid::unit(10, "pt"),
                                 ticks = TRUE)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 4)),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text  = element_text(size = 9),
      legend.background = element_rect(fill = scales::alpha("white", 0.9), color = NA),
      plot.margin = margin(2, 6, 2, 6)
    )
  
  # If absolutely all values are zero, optionally hide legend:
  if (vmax <= 0 || !is.finite(vmax)) p <- p + guides(fill = "none")
  p
}

# ── 5) Build panels (presence, not “only”) ────────────────────────
p_total <- make_map_from_vec(counts_presence %>% select(station, total_days),
                             "total_days", "Total Dust-Days",
                             option = "plasma",  trans = "sqrt")

p_du    <- make_map_from_vec(counts_presence %>% select(station, DU_days),
                             "DU_days",    "DU",
                             option = "viridis", trans = "sqrt")

p_bldu  <- make_map_from_vec(counts_presence %>% select(station, BLDU_days),
                             "BLDU_days",  "BLDU",
                             option = "magma",   trans = "sqrt")

p_ds    <- make_map_from_vec(counts_presence %>% select(station, DS_days),
                             "DS_days",    "DS",
                             option = "inferno", trans = "sqrt")

p_ss    <- make_map_from_vec(counts_presence %>% select(station, SS_days),
                             "SS_days",    "SS",
                             option = "turbo",   trans = "sqrt")

p_hz    <- make_map_from_vec(counts_presence %>% select(station, HZ_days),
                             "HZ_days",    "Dusty-HZ",
                             option = "cividis", trans = "sqrt")

p_multi <- make_map_from_vec(counts_presence %>% select(station, MULTI_PWC_days),
                             "MULTI_PWC_days", "Multiple PWC (≥2 codes/day)",
                             option = "rocket",  trans = "sqrt")

# Choose your layout (6 panels or include multi as 7th)
final_plot <- (p_total | p_du | p_bldu) /
  (p_ds    | p_ss | p_hz)

# If you want the Multi-PWC panel too, use:
# final_plot <- (p_total | p_du | p_bldu | p_ds) /
#               (p_ss    | p_hz | p_multi | plot_spacer())

final_plot
ggsave("dust_maps_presence.png", final_plot, width = 22, height = 10, dpi = 400)




#-----------------------------------------------------------
# ── Packages ─────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(viridis)
  library(ggrepel)
  library(maps)
  library(patchwork)
  library(lubridate)
  library(grid)   # unit()
})

# ── Inputs ───────────────────────────────────────────────────────
# Expect: dust_events_cv5 with columns: station, local_time (POSIXct), wxcodes, metar
CV_SHP <- "central_valley_alluvial_boundary_usgs.shp"

# Stations to plot (north→south list you’re using)
stations_df <- tibble::tribble(
  ~station, ~lat,      ~lon,
  "BFL",    35.43440, -119.0542,
  "BAB",    39.13609, -121.4366,
  "FAT",    36.78000, -119.7194,
  "HJO",    36.31139, -119.6232,
  "MAE",    36.98486, -120.1107,
  "MYV",    39.10203, -121.5688,
  "MCE",    37.28603, -120.5179,
  "OVE",    39.49000, -121.6200,
  "RBL",    40.15190, -122.2536,
  "RDD",    40.50900, -122.2934,
  "SAC",    38.50690, -121.4950,
  "SMF",    38.69542, -121.5908,
  "SCK",    37.89417, -121.2383,
  "VIS",    36.31867, -119.3929,
  "PTV",    36.02732, -119.0629
)

# ── Parse codes per report; collapse to daily presence ───────────
PWC_CODES <- c("BLDU","DU","DS","SS","HZ","VCBLDU")

norm_txt <- function(x) {
  x %>% stringr::str_to_upper() %>% stringr::str_replace_all("VC\\s*BLDU", "VCBLDU")
}
rx <- regex(paste0("\\b(", paste(PWC_CODES, collapse="|"), ")\\b"), ignore_case = TRUE)

per_report <- dust_events_cv5 %>%
  mutate(
    date    = as.Date(local_time),
    wxcodes = coalesce(as.character(wxcodes), ""),
    metar   = coalesce(as.character(metar),   ""),
    alltxt  = norm_txt(paste(wxcodes, metar, sep = " "))
  ) %>%
  mutate(
    codes_list  = str_extract_all(alltxt, rx) %>% map(~ unique(toupper(.x))),
    n_codes_row = map_int(codes_list, length),
    multi_row   = n_codes_row >= 2
  ) %>%
  select(station, date, codes_list, multi_row)

daily_presence <- per_report %>%
  unnest(codes_list, keep_empty = TRUE) %>%
  filter(!is.na(codes_list)) %>%
  mutate(code = codes_list) %>%
  distinct(station, date, code) %>%
  pivot_wider(names_from = code, values_from = code, values_fn = length, values_fill = 0)

for (cc in PWC_CODES) if (!cc %in% names(daily_presence)) daily_presence[[cc]] <- 0L

daily_presence <- daily_presence %>%
  transmute(
    station, date,
    DU_day     = DU      > 0,
    BLDU_day   = BLDU    > 0,
    VCBLDU_day = VCBLDU  > 0,
    DS_day     = DS      > 0,
    SS_day     = SS      > 0,
    HZ_day     = HZ      > 0
  ) %>%
  mutate(BLDU_day = BLDU_day & !VCBLDU_day)       # “plain” BLDU

# “≥2 distinct codes anywhere in the same day”
multi_day_anywhere <- per_report %>%
  mutate(code_set = map(codes_list, ~ setdiff(.x, character(0)))) %>%
  group_by(station, date) %>%
  summarise(n_distinct_codes = length(unique(unlist(code_set))), .groups = "drop") %>%
  mutate(MULTI_day = n_distinct_codes >= 2)

daily_flags <- stations_df %>% select(station) %>%
  left_join(daily_presence,     by = "station") %>%
  left_join(multi_day_anywhere, by = c("station","date")) %>%
  mutate(
    across(c(DU_day,BLDU_day,VCBLDU_day,DS_day,SS_day,HZ_day), ~ replace_na(.x, FALSE)),
    MULTI_day     = replace_na(MULTI_day, FALSE),
    ANY_dust_day  = DU_day | BLDU_day | VCBLDU_day | DS_day | SS_day | HZ_day
  )

# Per-station counts (presence)
counts_presence <- daily_flags %>%
  group_by(station) %>%
  summarise(
    total_days     = sum(ANY_dust_day, na.rm = TRUE),   # “dust days” per station
    DU_days        = sum(DU_day,       na.rm = TRUE),   # station-days with DU present
    BLDU_days      = sum(BLDU_day,     na.rm = TRUE),
    DS_days        = sum(DS_day,       na.rm = TRUE),
    SS_days        = sum(SS_day,       na.rm = TRUE),
    HZ_days        = sum(HZ_day,       na.rm = TRUE),
    MULTI_PWC_days = sum(MULTI_day,    na.rm = TRUE),   # station-days with ≥2 codes
    .groups = "drop"
  )

# ── Base maps ────────────────────────────────────────────────────
cv_boundary <- st_read(CV_SHP, quiet = TRUE) %>%
  st_zm(drop = TRUE) %>% st_transform(4326)

central_valley_counties <- c(
  "shasta","tehama","glenn","butte","colusa",
  "sutter","yuba","yolo","sacramento","san joaquin",
  "stanislaus","merced","madera","fresno","kings",
  "tulare","kern"
)
ca_counties <- map_data("county", "california") %>%
  filter(subregion %in% central_valley_counties)

# ── Map helper with explicit legend titles ───────────────────────
make_map_from_vec <- function(vals_df, value_col, title_key, legend_name,
                              option = "plasma", trans = "sqrt",
                              barheight_pt = 160) {
  
  dat <- stations_df %>%
    left_join(vals_df, by = "station") %>%
    mutate(value = !!rlang::sym(value_col), value = replace_na(value, 0L))
  
  vmax <- max(dat$value, na.rm = TRUE)
  vmax_adj <- if (is.finite(vmax) && vmax > 0) vmax else 1
  
  pos_df  <- dplyr::filter(dat, value > 0)
  zero_df <- dplyr::filter(dat, value == 0)
  sf_pos  <- sf::st_as_sf(pos_df,  coords = c("lon","lat"), crs = 4326)
  sf_zero <- sf::st_as_sf(zero_df, coords = c("lon","lat"), crs = 4326)
  
  p <- ggplot() +
    geom_sf(data = cv_boundary, fill = "#f0f0f0", color = "black", linewidth = 0.5) +
    geom_polygon(data = ca_counties, aes(long, lat, group = group),
                 fill = NA, color = "white", linewidth = 0.3) +
    geom_sf(data = sf_zero, shape = 21, fill = "white", color = "grey30", size = 6, linewidth = 0.6) +
    geom_sf(data = sf_pos,  aes(fill = value), shape = 21, color = "black", alpha = 0.9, size = 6) +
    scale_fill_viridis_c(
      name   = legend_name,
      option = option, begin = 0.15, end = 0.95,
      limits = c(0, vmax_adj), oob = scales::squish, trans = trans
    ) +
    geom_text_repel(data = dat, aes(lon, lat, label = station), size = 3, segment.size = 0.2) +
    coord_sf(
      xlim = sf::st_bbox(cv_boundary)[c("xmin","xmax")],
      ylim = sf::st_bbox(cv_boundary)[c("ymin","ymax")],
      expand = FALSE
    ) +
    labs(title = title_key) +
    theme_void(base_size = 12) +
    guides(fill = guide_colorbar(title.position = "top",
                                 direction = "vertical",
                                 barheight = grid::unit(barheight_pt, "pt"),
                                 barwidth  = grid::unit(10, "pt"),
                                 ticks = TRUE)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 4)),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text  = element_text(size = 9),
      legend.background = element_rect(fill = scales::alpha("white", 0.9), color = NA),
      plot.margin = margin(2, 6, 2, 6)
    )
  
  if (vmax <= 0 || !is.finite(vmax)) p <- p + guides(fill = "none")
  p
}

# ── Panels ───────────────────────────────────────────────────────
p_total <- make_map_from_vec(counts_presence %>% select(station, total_days),
                             "total_days", "Total Dust-Days",
                             legend_name = "days",
                             option = "plasma",  trans = "sqrt")

p_du    <- make_map_from_vec(counts_presence %>% select(station, DU_days),
                             "DU_days", "DU",
                             legend_name = "days",
                             option = "viridis", trans = "sqrt")

p_bldu  <- make_map_from_vec(counts_presence %>% select(station, BLDU_days),
                             "BLDU_days", "BLDU",
                             legend_name = "days",
                             option = "magma",   trans = "sqrt")

p_ds    <- make_map_from_vec(counts_presence %>% select(station, DS_days),
                             "DS_days", "DS",
                             legend_name = "days",
                             option = "inferno", trans = "sqrt")

p_ss    <- make_map_from_vec(counts_presence %>% select(station, SS_days),
                             "SS_days", "SS",
                             legend_name = "days",
                             option = "turbo",   trans = "sqrt")

p_hz    <- make_map_from_vec(counts_presence %>% select(station, HZ_days),
                             "HZ_days", "Dusty-HZ",
                             legend_name = "days",
                             option = "cividis", trans = "sqrt")

# Optional extra panel
p_multi <- make_map_from_vec(counts_presence %>% select(station, MULTI_PWC_days),
                             "MULTI_PWC_days", "Multiple PWC (≥2 codes/day)",
                             legend_name = "Station-days (≥2 codes)",
                             option = "rocket",  trans = "sqrt")

# Layout (6 panels; add p_multi if you want 7th/8th)
final_plot <- (p_total | p_du | p_bldu) /
  (p_ds    | p_ss | p_hz)

final_plot
ggsave("dust_maps_presence_legends.png", final_plot, width = 22, height = 10, dpi = 400)



#--------------------barplot2
# Packages
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(scales)

# -------------------- Build per-report flags --------------------
PWC_CODES <- c("BLDU","DU","DS","SS","HZ","VCBLDU")
norm_txt <- function(x) {
  x %>% stringr::str_to_upper() %>% stringr::str_replace_all("VC\\s*BLDU", "VCBLDU")
}
rx <- regex(paste0("\\b(", paste(PWC_CODES, collapse="|"), ")\\b"), ignore_case = TRUE)

per_report <- dust_events_cv5 %>%
  mutate(
    wxcodes = coalesce(as.character(wxcodes), ""),
    metar   = coalesce(as.character(metar),   ""),
    alltxt  = norm_txt(paste(wxcodes, metar, sep = " "))
  ) %>%
  mutate(
    codes_list  = str_extract_all(alltxt, rx) %>% map(~ unique(toupper(.x))),
    n_codes_row = map_int(codes_list, length),
    multi_row   = n_codes_row >= 2   # Multiple PWC in the SAME report
  ) %>%
  select(station, local_time, codes_list, multi_row)

report_flags <- per_report %>%
  mutate(
    has_DU     = map_lgl(codes_list, ~ "DU"     %in% .x),
    has_BLDU   = map_lgl(codes_list, ~ "BLDU"   %in% .x),
    has_VCBLDU = map_lgl(codes_list, ~ "VCBLDU" %in% .x),
    has_DS     = map_lgl(codes_list, ~ "DS"     %in% .x),
    has_SS     = map_lgl(codes_list, ~ "SS"     %in% .x),
    has_HZ     = map_lgl(codes_list, ~ "HZ"     %in% .x),
    BLDU_plain = has_BLDU & !has_VCBLDU,
    ANY_report = has_DU | BLDU_plain | has_DS | has_SS | has_HZ | has_VCBLDU
  )

# -------------------- Totals & percentages ----------------------
denom <- sum(report_flags$ANY_report, na.rm = TRUE)  # total dust reports
safe_pct <- function(n, d) if (is.finite(d) && d > 0) 100 * n / d else 0

cat_counts <- summarise(report_flags,
                        DU            = sum(has_DU,     na.rm = TRUE),
                        BLDU          = sum(BLDU_plain, na.rm = TRUE),
                        DS            = sum(has_DS,     na.rm = TRUE),
                        SS            = sum(has_SS,     na.rm = TRUE),
                        `Multiple PWC`= sum(multi_row,  na.rm = TRUE),
                        `Dusty-HZ`    = sum(has_HZ,     na.rm = TRUE),
                        .groups = "drop"
) %>%
  pivot_longer(everything(), names_to = "category", values_to = "count") %>%
  mutate(
    pct = safe_pct(count, denom),
    category = factor(category, levels = c("DU","BLDU","DS","SS","Multiple PWC","Dusty-HZ"))
  )

# -------------------- Plot (percent of dust reports) ------------
ggplot(cat_counts, aes(category, pct)) +
  geom_col(fill = "grey40", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            vjust = -0.5, size = 5) +
  scale_y_continuous(labels = percent_format(scale = 1),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL,
       y = "Share of dust reports (%)",
       title = "Composition of Dust Reports by Code (Report-level)") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(face = "bold")
  )

#-----------------------------------------------
# If needed:
if (!requireNamespace("ggpattern", quietly = TRUE)) install.packages("ggpattern")

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(scales)
library(ggpattern)

# ---------- Build per-report flags (same as before) ----------
PWC_CODES <- c("BLDU","DU","DS","SS","HZ","VCBLDU")
norm_txt <- function(x) x |>
  stringr::str_to_upper() |>
  stringr::str_replace_all("VC\\s*BLDU", "VCBLDU")
rx <- regex(paste0("\\b(", paste(PWC_CODES, collapse="|"), ")\\b"), ignore_case = TRUE)

per_report <- dust_events_cv5 %>%
  mutate(
    wxcodes = coalesce(as.character(wxcodes), ""),
    metar   = coalesce(as.character(metar),   ""),
    alltxt  = norm_txt(paste(wxcodes, metar, sep = " "))
  ) %>%
  mutate(
    codes_list  = str_extract_all(alltxt, rx) |> map(~ unique(toupper(.x))),
    n_codes_row = map_int(codes_list, length),
    multi_row   = n_codes_row >= 2
  ) %>%
  select(station, local_time, codes_list, multi_row)

report_flags <- per_report %>%
  mutate(
    has_DU     = map_lgl(codes_list, ~ "DU"     %in% .x),
    has_BLDU   = map_lgl(codes_list, ~ "BLDU"   %in% .x),
    has_VCBLDU = map_lgl(codes_list, ~ "VCBLDU" %in% .x),
    has_DS     = map_lgl(codes_list, ~ "DS"     %in% .x),
    has_SS     = map_lgl(codes_list, ~ "SS"     %in% .x),
    has_HZ     = map_lgl(codes_list, ~ "HZ"     %in% .x),
    BLDU_plain = has_BLDU & !has_VCBLDU,  # typo fixed below line!
    BLDU_plain = has_BLDU & !has_VCBLDU,
    ANY_report = has_DU | BLDU_plain | has_DS | has_SS | has_HZ | has_VCBLDU
  )

denom <- sum(report_flags$ANY_report, na.rm = TRUE)
safe_pct <- function(n, d) if (is.finite(d) && d > 0) 100 * n / d else 0

cat_counts <- summarise(report_flags,
                        DU             = sum(has_DU,     na.rm = TRUE),
                        BLDU           = sum(BLDU_plain, na.rm = TRUE),
                        DS             = sum(has_DS,     na.rm = TRUE),
                        SS             = sum(has_SS,     na.rm = TRUE),
                        `Multiple PWC` = sum(multi_row,  na.rm = TRUE),
                        `Dusty-HZ`     = sum(has_HZ,     na.rm = TRUE),
                        .groups = "drop"
) %>%
  pivot_longer(everything(), names_to = "category", values_to = "count") %>%
  mutate(
    pct = safe_pct(count, denom),
    category = factor(category, levels = c("DU","BLDU","DS","SS","Multiple PWC","Dusty-HZ"))
  )

# ---------- Colors for fill and for the stripe pattern ----------
pal_fill <- c(
  "DU"            = "#4C78A8",
  "BLDU"          = "#B279A2",
  "DS"            = "#F58518",
  "SS"            = "#54A24B",
  "Multiple PWC"  = "#E45756",
  "Dusty-HZ"      = "#808080"
)
pal_stripe <- c(
  "DU"            = "#1B3A6B",
  "BLDU"          = "#6B3F7A",
  "DS"            = "#9D3F0E",
  "SS"            = "#2E6F31",
  "Multiple PWC"  = "#8C1D18",
  "Dusty-HZ"      = "#4A4A4A"
)

cat_counts <- cat_counts %>%
  mutate(fill_col = pal_fill[category],
         stripe_col = pal_stripe[category])

# ---------- Plot with ggpattern (/// stripes) ----------
ggplot(cat_counts, aes(category, pct)) +
  geom_col_pattern(
    aes(fill = category,
        pattern_fill   = stripe_col,
        pattern_colour = stripe_col),
    pattern          = "stripe",
    pattern_angle    = 60,     # “///”
    pattern_spacing  = 0.05,   # adjust to taste (smaller = denser)
    pattern_density  = 0.5,
    width            = 0.7,
    colour           = "black"
  ) +
  scale_fill_manual(values = pal_fill, guide = "none") +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 5) +
  scale_y_continuous(labels = percent_format(scale = 1),
                     expand  = expansion(mult = c(0, 0.1))) +
  labs(x = NULL,
       y = "Share of dust reports (%)",
       title = "Composition of Dust Reports by Code (Report-level, striped)") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(face = "bold")
  )

#-----------------------------------------
# ── Packages ─────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(ggplot2); library(scales)
})
if (!requireNamespace("ggpattern", quietly = TRUE)) install.packages("ggpattern")
library(ggpattern)

# INPUT: dust_events_cv5 with columns: station, local_time (POSIXct), wxcodes, metar

# ── 1) Extract codes per report → station/date/code presence ─────
PWC <- c("BLDU","VCBLDU","DU","DS","SS","HZ")
norm_txt <- function(x) x |> str_to_upper() |> str_replace_all("VC\\s*BLDU","VCBLDU")
rx <- paste0("\\b(", paste(PWC, collapse="|"), ")\\b")

daily_long <- dust_events_cv5 %>%
  transmute(
    station,
    date = as.Date(local_time),
    txt  = norm_txt(paste(coalesce(as.character(wxcodes), ""),
                          coalesce(as.character(metar),   "")))
  ) %>%
  mutate(code = str_extract_all(txt, rx)) %>%
  unnest_longer(code, values_to = "code") %>%
  filter(!is.na(code)) %>%
  mutate(code = toupper(code)) %>%
  distinct(station, date, code)                       # presence, not counts

# ── 2) Build day flags; BLDU excludes VCBLDU; MULTI = ≥2 distinct codes/day ─
daily_wide <- daily_long %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = code, values_from = present, values_fill = 0L)

for (cc in PWC) if (!cc %in% names(daily_wide)) daily_wide[[cc]] <- 0L

daily_flags <- daily_wide %>%
  mutate(across(all_of(PWC), ~ . > 0)) %>%
  mutate(
    BLDU_plain  = BLDU & !VCBLDU,
    DU_day      = DU,
    BLDU_day    = BLDU_plain,
    DS_day      = DS,
    SS_day      = SS,
    HZ_day      = HZ,
    VCBLDU_day  = VCBLDU,
    MULTI_day   = (DU_day + BLDU_day + VCBLDU_day + DS_day + SS_day + HZ_day) >= 2
  ) %>%
  select(station, date, DU_day, BLDU_day, DS_day, SS_day, HZ_day, MULTI_day)

# ── 3) Totals (day-level counts) and % of TOTAL counts across categories ─
cat_counts_day <- tibble::tibble(
  category = factor(c("DU","BLDU","DS","SS","Multiple PWC","Dusty-HZ"),
                    levels = c("DU","BLDU","DS","SS","Multiple PWC","Dusty-HZ")),
  count = c(
    sum(daily_flags$DU_day,    na.rm = TRUE),
    sum(daily_flags$BLDU_day,  na.rm = TRUE),
    sum(daily_flags$DS_day,    na.rm = TRUE),
    sum(daily_flags$SS_day,    na.rm = TRUE),
    sum(daily_flags$MULTI_day, na.rm = TRUE),
    sum(daily_flags$HZ_day,    na.rm = TRUE)
  )
)

cat_counts_total <- cat_counts_day %>%
  mutate(
    pct_of_total     = count / sum(count),
    pct_of_total_lab = percent(pct_of_total, accuracy = 0.1)
  )

print(cat_counts_total)

# ── 4) Plot: % of total counts (striped bars, one legend "Code") ─
pal_fill <- c("DU"="#4C78A8","BLDU"="#B279A2","DS"="#F58518","SS"="#54A24B",
              "Multiple PWC"="#E45756","Dusty-HZ"="#808080")
pal_stripe <- c("DU"="#1B3A6B","BLDU"="#6B3F7A","DS"="#9D3F0E","SS"="#2E6F31",
                "Multiple PWC"="#8C1D18","Dusty-HZ"="#4A4A4A")

ggplot(cat_counts_total, aes(category, pct_of_total, fill = category)) +
  ggpattern::geom_col_pattern(
    aes(pattern_fill = category, pattern_colour = category),
    pattern = "stripe", pattern_angle = 60, pattern_spacing = 0.05,
    pattern_density = 0.5, width = 0.7, colour = "black"
  ) +
  scale_fill_manual(name = "Code", values = pal_fill) +
  scale_pattern_fill_manual(values = pal_stripe, guide = "none") +
  scale_pattern_colour_manual(values = pal_stripe, guide = "none") +
  guides(fill = guide_legend(override.aes = list(pattern = "stripe"))) +
  geom_text(aes(label = pct_of_total_lab), vjust = -0.5, size = 5) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Percentage of total Dust codes (%)",
       title = "") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(face = "bold"),
        legend.position = "right")



#---------------------------year to yar variability of dust code plot below
# ── Packages ─────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(lubridate); library(ggplot2); library(scales)
})
if (!requireNamespace("ggpattern", quietly = TRUE)) install.packages("ggpattern")
library(ggpattern)

# INPUT: dust_events_cv5 with columns: station, local_time (POSIXct), wxcodes, metar

# ── 1) Extract codes per report → station/date/code presence ─────
# Codes tracked; we’ll treat VCBLDU as BLDU family for categorization
PWC <- c("BLDU","VCBLDU","DU","DS","SS","HZ")
norm_txt <- function(x) x |>
  str_to_upper() |>
  str_replace_all("VC\\s*BLDU","VCBLDU")
rx <- paste0("\\b(", paste(PWC, collapse="|"), ")\\b")

daily_long <- dust_events_cv5 %>%
  transmute(
    station,
    date = as.Date(local_time),
    txt  = norm_txt(paste(coalesce(as.character(wxcodes), ""),
                          coalesce(as.character(metar),   "")))
  ) %>%
  mutate(code = str_extract_all(txt, rx)) %>%
  unnest_longer(code, values_to = "code") %>%
  filter(!is.na(code)) %>%
  mutate(code = toupper(code)) %>%
  distinct(station, date, code)        # presence (not hourly counts)

# ── 2) Build station–day flags; BLDU family; Multiple PWC ────────
daily_wide <- daily_long %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = code, values_from = present, values_fill = 0L)

# make sure all exist
for (cc in PWC) if (!cc %in% names(daily_wide)) daily_wide[[cc]] <- 0L

daily_flags <- daily_wide %>%
  mutate(across(all_of(PWC), ~ . > 0)) %>%
  mutate(
    DU_day      = DU,
    BLDUfam_day = BLDU | VCBLDU,         # BLDU family (includes VCBLDU for this figure)
    DS_day      = DS,
    SS_day      = SS,
    HZ_day      = HZ
  ) %>%
  # how many families active today (BLDU + VCBLDU count as one family)
  mutate(
    n_families  = DU_day + BLDUfam_day + DS_day + SS_day + HZ_day,
    MULTI_day   = n_families >= 2,
    ANY_day     = n_families >= 1
  ) %>%
  select(station, date, DU_day, BLDUfam_day, DS_day, SS_day, HZ_day, MULTI_day, ANY_day)

# ── 3) Assign a single, exclusive category per station–day ───────
# Priority: MULTI → DU → BLDU → DS → SS → Dusty-HZ (HZ only)
daily_cat <- daily_flags %>%
  mutate(
    category = case_when(
      MULTI_day                 ~ "Multiple PWC",
      DU_day                    ~ "DU",
      BLDUfam_day               ~ "BLDU",
      DS_day                    ~ "DS",
      SS_day                    ~ "SS",
      HZ_day                    ~ "Dusty-HZ",
      TRUE                      ~ NA_character_
    )
  ) %>%
  filter(!is.na(category)) %>%
  mutate(year = year(date)) %>%
  filter(year >= 2005, year <= 2024)

# ── 4) Yearly counts and in-year % fractions (stack to 100%) ─────
yearly_frac <- daily_cat %>%
  count(year, category, name = "count") %>%
  group_by(year) %>%
  mutate(pct = 100 * count / sum(count)) %>%
  ungroup() %>%
  # ensure all categories exist each year (fill 0)
  complete(
    year,
    category = factor(c("DU","BLDU","DS","SS","Multiple PWC","Dusty-HZ"),
                      levels = c("Dusty-HZ","SS","DS","BLDU","DU","Multiple PWC")), # control stack order
    fill = list(count = 0, pct = 0)
  ) %>%
  arrange(year, category)

# ── 5) Plot: 100% stacked bars (with stripe on Multiple PWC) ─────
pal_fill <- c(
  "DU"            = "#A86B00",
  "BLDU"          = "#7DB1FF",
  "DS"            = "#1B7F6A",
  "SS"            = "#52B788",
  "Multiple PWC"  = "#444444",
  "Dusty-HZ"      = "#EAD97F"
)

ggplot(yearly_frac,
       aes(x = factor(year), y = pct, fill = category)) +
  ggpattern::geom_col_pattern(
    aes(pattern = ifelse(category == "Multiple PWC", "stripe", "none"),
        pattern_angle = ifelse(category == "Multiple PWC", 45, 0),
        pattern_spacing = ifelse(category == "Multiple PWC", 0.04, 0.02),
        pattern_density = ifelse(category == "Multiple PWC", 0.5, 0)),
    width = 0.9, color = "black", linewidth = 0.2
  ) +
  scale_fill_manual(values = pal_fill, name = "Dust Category") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, 100, 20)) +
  labs(x = NULL, y = "Percentage Contribution (%)",
       title = "Percentage fraction of identified Dust Category (2005–2024)") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )


# ── Packages ─────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(lubridate); library(ggplot2); library(scales)
})
if (!requireNamespace("ggpattern", quietly = TRUE)) install.packages("ggpattern")
library(ggpattern)

# INPUT expected in memory: dust_events_cv5 with columns:
#   station, local_time (POSIXct), wxcodes, metar

# ── 1) Per-report codes → station–day presence ───────────────────
PWC <- c("BLDU","VCBLDU","DU","DS","SS","HZ")

norm_txt <- function(x) {
  x |>
    stringr::str_to_upper() |>
    stringr::str_replace_all("VC\\s*BLDU", "VCBLDU")
}
rx <- paste0("\\b(", paste(PWC, collapse="|"), ")\\b")

daily_long <- dust_events_cv5 %>%
  transmute(
    station,
    date = as.Date(local_time),
    txt  = norm_txt(paste(coalesce(as.character(wxcodes), ""),
                          coalesce(as.character(metar),   "")))
  ) %>%
  mutate(code = str_extract_all(txt, rx)) %>%
  unnest_longer(code, values_to = "code") %>%
  filter(!is.na(code)) %>%
  mutate(code = toupper(code)) %>%
  distinct(station, date, code)     # presence (not hourly counts)

# ── 2) Station–day flags; BLDU family; Multiple PWC ──────────────
daily_wide <- daily_long %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = code, values_from = present, values_fill = 0L)

# Ensure all needed columns exist
for (cc in PWC) if (!cc %in% names(daily_wide)) daily_wide[[cc]] <- 0L

daily_flags <- daily_wide %>%
  mutate(across(all_of(PWC), ~ . > 0)) %>%
  mutate(
    DU_day      = DU,
    BLDUfam_day = BLDU | VCBLDU,          # treat VCBLDU as BLDU family (for this figure)
    DS_day      = DS,
    SS_day      = SS,
    HZ_day      = HZ,
    n_families  = DU_day + BLDUfam_day + DS_day + SS_day + HZ_day,
    MULTI_day   = n_families >= 2,
    ANY_day     = n_families >= 1
  ) %>%
  select(station, date, DU_day, BLDUfam_day, DS_day, SS_day, HZ_day, MULTI_day, ANY_day)

# ── 3) Assign ONE exclusive category per station–day ─────────────
# Priority: MULTI → DU → BLDU → DS → SS → Dusty-HZ (HZ only)
daily_cat <- daily_flags %>%
  mutate(
    category = case_when(
      MULTI_day     ~ "Multiple PWC",
      DU_day        ~ "DU",
      BLDUfam_day   ~ "BLDU",
      DS_day        ~ "DS",
      SS_day        ~ "SS",
      HZ_day        ~ "Dusty-HZ",
      TRUE          ~ NA_character_
    ),
    year = lubridate::year(date)
  ) %>%
  filter(!is.na(category), year >= 2005, year <= 2024)

# ── 4) Yearly counts and in-year % fractions (sum to 100%) ───────
cat_levels <- c("Dusty-HZ","SS","DS","BLDU","DU","Multiple PWC")  # stack order bottom→top

yearly_frac <- daily_cat %>%
  count(year, category, name = "count") %>%
  group_by(year) %>%
  mutate(pct = 100 * count / sum(count)) %>%
  ungroup() %>%
  complete(year,
           category = factor(cat_levels, levels = cat_levels),
           fill = list(count = 0, pct = 0)) %>%
  arrange(year, category)

# ── 5) Plot: 100% stacked with ONE legend (categories only) ──────
pal_fill <- c(
  "DU"            = "#A86B00",
  "BLDU"          = "#7DB1FF",
  "DS"            = "#1B7F6A",
  "SS"            = "#52B788",
  "Multiple PWC"  = "#444444",
  "Dusty-HZ"      = "#EAD97F"
)

p <- ggplot(yearly_frac,
            aes(x = factor(year), y = pct, fill = category)) +
  ggpattern::geom_col_pattern(
    # map ONLY pattern on/off to avoid extra legend keys
    aes(pattern = ifelse(category == "Multiple PWC", "stripe", "none")),
    pattern_angle   = 45,
    pattern_spacing = 0.04,
    pattern_density = 0.5,
    width = 0.9, color = "black", linewidth = 0.2
  ) +
  scale_fill_manual(values = pal_fill, name = "Dust Category") +
  # hide pattern guides completely
  ggpattern::scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"),
                                  guide = "none") +
  guides(pattern = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Percentage Contribution (%)",
       title = "Percentage fraction of identified Dust Category (2005–2024)") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

p
ggsave("dust_category_share_by_year.png", p, width = 14, height = 7, dpi = 300)

#--------------------------------------
# ── 5) Plot: 100% stacked with one legend (colors + slashes) ─────

# Color palette (clear separation for DS vs SS)
pal_fill <- c(
  "DU"            = "#A86B00",  # brown/orange
  "BLDU"          = "#7DB1FF",  # light blue
  "DS"            = "#009E73",  # green (Okabe–Ito)
  "SS"            = "#CC79A7",  # magenta (Okabe–Ito)  << more distinct from DS
  "Multiple PWC"  = "#585858",  # dark grey
  "Dusty-HZ"      = "#EAD97F"   # pale gold
)

# Pattern map: only "Multiple PWC" is striped
pat_vals <- c(
  "DU"            = "none",
  "BLDU"          = "none",
  "DS"            = "none",
  "SS"            = "none",
  "Multiple PWC"  = "stripe",
  "Dusty-HZ"      = "none"
)

p <- ggplot(yearly_frac,
            aes(x = factor(year),
                y = pct,
                fill = category,
                pattern = category)) +   # map pattern to category
  ggpattern::geom_col_pattern(
    width = 0.9, color = "black", linewidth = 0.2,
    pattern_angle   = 45,
    pattern_spacing = 0.04,
    pattern_density = 0.5,
    pattern_fill    = "grey20",
    pattern_colour  = "grey20"
  ) +
  # Give BOTH scales the same legend name so they merge into one
  scale_fill_manual(values = pal_fill, name = "Dust Category") +
  ggpattern::scale_pattern_manual(values = pat_vals, name = "Dust Category") +
  guides(
    # slightly thicker key border helps show the pattern
    fill    = guide_legend(order = 1, override.aes = list(colour = "black", linewidth = 0.3)),
    pattern = guide_legend(order = 1, override.aes = list(colour = "black", linewidth = 0.3))
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Percentage Contribution (%)",
       title = "Percentage fraction of identified Dust Category (2005–2024)") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 90, vjust = 0.5),
    legend.position    = "top",
    legend.title       = element_text(face = "bold")
  )

p
ggsave("dust_category_share_by_year.png", p, width = 14, height = 7, dpi = 300)




# ── 5) Plot: 100% stacked with one legend (colors + slashes) ─────

# Color palette (clear separation for DS vs SS)
pal_fill <- c(
  "DU"            = "#A86B00",  # brown/orange
  "BLDU"          = "#7DB1FF",  # light blue
  "DS"            = "#009E73",  # green (Okabe–Ito)
  "SS"            = "#CC79A7",  # magenta (Okabe–Ito)  << more distinct from DS
  "Multiple PWC"  = "#585858",  # dark grey
  "Dusty-HZ"      = "#EAD97F"   # pale gold
)

# Pattern map: only "Multiple PWC" is striped
pat_vals <- c(
  "DU"            = "none",
  "BLDU"          = "none",
  "DS"            = "none",
  "SS"            = "none",
  "Multiple PWC"  = "stripe",
  "Dusty-HZ"      = "none"
)

p <- ggplot(yearly_frac,
            aes(x = factor(year),
                y = pct,
                fill = category,
                pattern = category)) +   # map pattern to category
  ggpattern::geom_col_pattern(
    width = 0.9, color = "black", linewidth = 0.2,
    pattern_angle   = 45,
    pattern_spacing = 0.04,
    pattern_density = 0.5,
    pattern_fill    = "grey20",
    pattern_colour  = "grey20"
  ) +
  # Give BOTH scales the same legend name so they merge into one
  scale_fill_manual(values = pal_fill, name = "Dust Category") +
  ggpattern::scale_pattern_manual(values = pat_vals, name = "Dust Category") +
  guides(
    # slightly thicker key border helps show the pattern
    fill    = guide_legend(order = 1, override.aes = list(colour = "black", linewidth = 0.3)),
    pattern = guide_legend(order = 1, override.aes = list(colour = "black", linewidth = 0.3))
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Percentage fraction (%) of identified Dust Category",
       title = "") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 90, vjust = 0.5),
    legend.position    = "top",
    legend.title       = element_text(face = "bold")
  )

p
ggsave("share_by_year.png", p, width = 14, height = 7, dpi = 300)


#---------------KEY----------------------------------------------------------------------


# ── Packages ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(lubridate)
  library(ggplot2); library(scales); library(patchwork)
})
if (!requireNamespace("ggpattern", quietly = TRUE)) install.packages("ggpattern")
library(ggpattern)

# INPUT: dust_events_cv5 with columns: station, local_time (POSIXct), wxcodes, metar

# ── 1) Parse codes per report → station–day presence ─────────────────────
PWC <- c("BLDU","VCBLDU","DU","DS","SS","HZ")
norm_txt <- function(x) x |>
  stringr::str_to_upper() |>
  stringr::str_replace_all("VC\\s*BLDU", "VCBLDU")
rx <- paste0("\\b(", paste(PWC, collapse="|"), ")\\b")

daily_long <- dust_events_cv5 %>%
  transmute(
    station,
    date = as.Date(local_time),
    txt  = norm_txt(paste(coalesce(as.character(wxcodes), ""),
                          coalesce(as.character(metar),   "")))
  ) %>%
  mutate(code = stringr::str_extract_all(txt, rx)) %>%
  tidyr::unnest_longer(code) %>%
  filter(!is.na(code)) %>%
  mutate(code = toupper(code)) %>%
  distinct(station, date, code)  # presence (not hourly counts)

# Wide flags (ensure all columns exist)
daily_wide <- daily_long %>%
  mutate(present = 1L) %>%
  tidyr::pivot_wider(names_from = code, values_from = present, values_fill = 0L)
for (cc in PWC) if (!cc %in% names(daily_wide)) daily_wide[[cc]] <- 0L

# Family/day flags + Multi
daily_flags <- daily_wide %>%
  mutate(across(all_of(PWC), ~ . > 0)) %>%
  transmute(
    station, date,
    DU_day      = DU,
    BLDUfam_day = BLDU | VCBLDU,     # BLDU family
    DS_day      = DS,
    SS_day      = SS,
    HZ_day      = HZ,
    n_families  = (DU | (BLDU | VCBLDU) | DS | SS | HZ) + 0L, # bool->int then sum via "+"
    n_families  = (DU_day + BLDUfam_day + DS_day + SS_day + HZ_day),
    MULTI_day   = n_families >= 2,
    ANY_day     = n_families >= 1
  )

# ── 2) Exclusive category per station–day (priority) ────────────────────
# Priority: MULTI → DU → BLDU → DS → SS → Dusty-HZ
daily_cat <- daily_flags %>%
  mutate(
    category = case_when(
      MULTI_day     ~ "Multiple PWC",
      DU_day        ~ "DU",
      BLDUfam_day   ~ "BLDU",
      DS_day        ~ "DS",
      SS_day        ~ "SS",
      HZ_day        ~ "Dusty-HZ",
      TRUE          ~ NA_character_
    ),
    year = lubridate::year(date)
  ) %>%
  filter(!is.na(category), dplyr::between(year, 2005, 2024))

cat_levels <- c("Dusty-HZ","SS","DS","BLDU","DU","Multiple PWC") # bottom→top order
pal_fill <- c(
  "DU"            = "#A86B00",
  "BLDU"          = "#7DB1FF",
  "DS"            = "#1B7F6A",
  "SS"            = "#35C39B",
  "Multiple PWC"  = "#4A4A4A",
  "Dusty-HZ"      = "#EAD97F"
)

# ── 3) Yearly percentages (composition) ─────────────────────────────────
yearly_frac <- daily_cat %>%
  count(year, category, name = "count") %>%
  group_by(year) %>% mutate(pct = 100 * count / sum(count)) %>% ungroup() %>%
  complete(year = 2005:2024,
           category = factor(cat_levels, levels = cat_levels),
           fill = list(count = 0, pct = 0)) %>%
  mutate(year = factor(year, levels = 2005:2024)) %>%
  arrange(year, category)

# ── 4) Yearly absolute counts (+ totals) ────────────────────────────────
yearly_counts <- daily_cat %>%
  count(year, category, name = "count") %>%
  complete(year = 2005:2024,
           category = factor(cat_levels, levels = cat_levels),
           fill = list(count = 0)) %>%
  mutate(year = factor(year, levels = 2005:2024)) %>%
  arrange(year, category)

totals_per_year <- yearly_counts %>%
  group_by(year) %>% summarise(total = sum(count), .groups = "drop")

ymax_counts <- max(totals_per_year$total, na.rm = TRUE)
ybreaks_cnt  <- pretty(c(0, ymax_counts))

# ── 5) Plots ────────────────────────────────────────────────────────────
# (A) Composition (100%)
p_comp <- ggplot(yearly_frac,
                 aes(x = year, y = pct, fill = category)) +
  ggpattern::geom_col_pattern(
    aes(pattern = ifelse(category == "Multiple PWC", "stripe", "none")),
    pattern_angle   = 45,
    pattern_spacing = 0.04,
    pattern_density = 0.5,
    width = 0.9, color = "black", linewidth = 0.2,
    show.legend = TRUE
  ) +
  scale_fill_manual(values = pal_fill, name = "Dust Category") +
  ggpattern::scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"),
                                  guide = "none") +
  scale_y_continuous(labels = label_percent(scale = 1),
                     breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Percentage Contribution (%)",
       title = "(a) % Fraction/year ") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

# (B) Absolute counts (stacked) + total line
p_abs <- ggplot(yearly_counts,
                aes(x = year, y = count, fill = category)) +
  ggpattern::geom_col_pattern(
    aes(pattern = ifelse(category == "Multiple PWC", "stripe", "none")),
    pattern_angle   = 45,
    pattern_spacing = 0.04,
    pattern_density = 0.5,
    width = 0.9, color = "black", linewidth = 0.2,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = pal_fill) +
  ggpattern::scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"),
                                  guide = "none") +
  geom_line(data = totals_per_year, aes(x = year, y = total, group = 1),
            color = "black", linewidth = 0.7, inherit.aes = FALSE) +
  geom_point(data = totals_per_year, aes(x = year, y = total),
             color = "black", size = 1.6, inherit.aes = FALSE) +
  scale_y_continuous(breaks = ybreaks_cnt, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Dusty-days",
       title = "(b) Annual Dusty-days") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

# Combine
paired_plot <- p_comp / p_abs + plot_layout(heights = c(1, 1.05))
paired_plot

ggsave("dust_category_composition_and_counts.png", paired_plot,
       width = 14, height = 11, dpi = 300)




#------------------plot composite-----------------------

# ── Packages ─────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(lubridate); library(purrr)
  library(sf); library(ggplot2); library(ggrepel)
  library(viridis); library(maps); library(patchwork)
})

# ── Inputs ───────────────────────────────────────────────────────
# Required in memory: dust_events_cv5 with at least
# station, local_time (POSIXct), wxcodes, metar, tmpc,
# relh (or RH_calc), wind_speed_mps, visibility_km (or vsby)
CV_SHP <- "central_valley_alluvial_boundary_usgs.shp"

# 15 stations (no FCH)
stations_df <- tibble::tribble(
  ~station, ~lat,      ~lon,
  "BFL",    35.43440, -119.0542,
  "BAB",    39.13609, -121.4366,
  "FAT",    36.78000, -119.7194,
  "HJO",    36.31139, -119.6232,
  "MAE",    36.98486, -120.1107,
  "MYV",    39.10203, -121.5688,
  "MCE",    37.28603, -120.5179,
  "OVE",    39.49000, -121.6200,
  "RBL",    40.15190, -122.2536,
  "RDD",    40.50900, -122.2934,
  "SAC",    38.50690, -121.4950,
  "SMF",    38.69542, -121.5908,
  "SCK",    37.89417, -121.2383,
  "VIS",    36.31867, -119.3929,
  "PTV",    36.02732, -119.0629
)

# ── Helpers ──────────────────────────────────────────────────────
PWC_CODES <- c("DU","BLDU","VCBLDU","DS","SS","HZ")
rx_codes   <- regex(paste0("\\b(", paste(PWC_CODES, collapse="|"), ")\\b"),
                    ignore_case = TRUE)

norm_txt <- function(x) {
  x %>% stringr::str_to_upper() %>% stringr::str_replace_all("VC\\s*BLDU", "VCBLDU")
}

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# Convenience flags for optional columns (so mutate won't touch missing ones)
has_relh <- "relh"          %in% names(dust_events_cv5)
has_rhc  <- "RH_calc"       %in% names(dust_events_cv5)
has_visk <- "visibility_km" %in% names(dust_events_cv5)
has_vsby <- "vsby"          %in% names(dust_events_cv5)

# ── 1) Identify dust days (any PWC code present that date) ───────
dust_days <- dust_events_cv5 %>%
  transmute(
    station,
    date = as.Date(local_time),
    txt  = norm_txt(paste(coalesce(as.character(wxcodes), ""),
                          coalesce(as.character(metar),   "")))
  ) %>%
  mutate(any_pwc = str_detect(txt, rx_codes)) %>%
  group_by(station, date) %>%
  summarise(dust_day = any(any_pwc, na.rm = TRUE), .groups = "drop") %>%
  filter(dust_day) %>%
  semi_join(stations_df, by = "station")

# ── 2) Daily means on dust days (per station–date) ───────────────
daily_means <- dust_events_cv5 %>%
  mutate(date = as.Date(local_time)) %>%
  inner_join(dust_days, by = c("station","date")) %>%
  mutate(
    temp_c = suppressWarnings(as.numeric(tmpc)),
    
    # RH: prefer relh, else RH_calc (using dplyr 1.1+ pick/any_of safely)
    rh_pct = if (has_relh || has_rhc) {
      do.call(coalesce, c(pick(any_of(c("relh","RH_calc"))), list(NA_real_)))
    } else NA_real_,
    
    wind_m = suppressWarnings(as.numeric(wind_speed_mps)),
    
    # Visibility: prefer visibility_km, else convert vsby (miles->km)
    vis_km = dplyr::case_when(
      has_visk             ~ suppressWarnings(as.numeric(.data$visibility_km)),
      !has_visk & has_vsby ~ suppressWarnings(as.numeric(.data$vsby) * 1.60934),
      TRUE                 ~ NA_real_
    )
  ) %>%
  group_by(station, date) %>%
  summarise(
    temp_c_day = safe_mean(temp_c),
    rh_day     = safe_mean(rh_pct),
    wind_day   = safe_mean(wind_m),
    vis_day    = safe_mean(vis_km),
    .groups = "drop"
  )

# ── 3) Composite (mean across dust days) per station ─────────────
comp_stats <- daily_means %>%
  group_by(station) %>%
  summarise(
    n_dust_days = n(),
    temp_c      = safe_mean(temp_c_day),
    rh_pct      = safe_mean(rh_day),
    wind_mps    = safe_mean(wind_day),
    vis_km      = safe_mean(vis_day),
    .groups = "drop"
  ) %>%
  right_join(stations_df, by = "station")   # keep all 15 stations (NA if no dust days)

# ── 4) Geometry: CV boundary & counties for context ──────────────
cv_boundary <- st_read(CV_SHP, quiet = TRUE) %>%
  st_zm(drop = TRUE) %>%
  st_transform(4326)

central_valley_counties <- c(
  "shasta","tehama","glenn","butte","colusa",
  "sutter","yuba","yolo","sacramento","san joaquin",
  "stanislaus","merced","madera","fresno","kings",
  "tulare","kern"
)
ca_counties <- map_data("county", "california") %>%
  filter(subregion %in% central_valley_counties)

pts_sf <- st_as_sf(comp_stats, coords = c("lon","lat"), crs = 4326)

# ── 5) Small panel factory ───────────────────────────────────────
panel_map <- function(sf_pts, value_col, title_str, legend_title,
                      palette = "magma", trans = "identity") {
  
  vals <- sf_pts[[value_col]]
  rng  <- range(vals, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    rng <- c(0, 1)                  # safe default if degenerate
  } else {
    pad <- 0.05 * diff(rng); rng <- c(rng[1] - pad, rng[2] + pad)
  }
  
  show_pts  <- is.finite(vals) & !is.na(vals)
  sf_pos    <- sf_pts[ show_pts, ]
  sf_hollow <- sf_pts[!show_pts, ]
  
  ggplot() +
    geom_sf(data = cv_boundary, fill = "#f0f0f0", color = "black", linewidth = 0.5) +
    geom_polygon(data = ca_counties, aes(long, lat, group = group),
                 fill = NA, color = "white", linewidth = 0.3) +
    geom_sf(data = sf_hollow, shape = 21, fill = "white", color = "grey30", size = 4) +
    geom_sf(data = sf_pos, aes(fill = .data[[value_col]]),
            shape = 21, color = "black", size = 4) +
    scale_fill_viridis_c(
      name   = legend_title,
      option = palette,
      limits = rng, oob = scales::squish,
      trans  = trans
    ) +
    geom_text_repel(
      data = st_drop_geometry(sf_pts),
      aes(x = st_coordinates(sf_pts)[,1], y = st_coordinates(sf_pts)[,2], label = station),
      size = 3, min.segment.length = 0.05
    ) +
    coord_sf(
      xlim = st_bbox(cv_boundary)[c("xmin","xmax")],
      ylim = st_bbox(cv_boundary)[c("ymin","ymax")],
      expand = FALSE
    ) +
    labs(title = title_str) +
    theme_void(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
}

# ── 6) Build the four panels ─────────────────────────────────────
p_temp <- panel_map(pts_sf, "temp_c",   "(a) Temperature (°C)", "Temp (°C)",       palette = "plasma")
p_rh   <- panel_map(pts_sf, "rh_pct",   "(b) RH (%)",           "RH (%)",          palette = "magma")
p_wind <- panel_map(pts_sf, "wind_mps", "(c) Wind speed (m/s)", "Wind (m/s)",      palette = "inferno")
p_vis  <- panel_map(pts_sf, "vis_km",   "(d) Visibility (km)",  "Visibility (km)", palette = "viridis")

final_plot <- (p_temp | p_rh | p_wind | p_vis)
final_plot

final_plot <- (p_rh | p_wind | p_vis)

# Optional: save figure
ggsave("composite_surface_on_dust_days.png", final_plot, width = 12, height = 10, dpi = 400)

# --- If you prefer fixed colorbar ranges, replace inside panel_map() per variable, e.g.:
# p_temp <- panel_map(pts_sf, "temp_c", "(a) Temperature (°C)", "Temp (°C)") + 
#           scale_fill_viridis_c(limits = c(18, 28), option = "plasma")
# p_rh   <- panel_map(pts_sf, "rh_pct", "(b) RH (%)", "RH (%)") + 
#           scale_fill_viridis_c(limits = c(25, 45), option = "magma")
# p_wind <- panel_map(pts_sf, "wind_mps", "(c) Wind speed (m/s)", "Wind (m/s)") + 
#           scale_fill_viridis_c(limits = c(7, 11), option = "inferno")
# p_vis  <- panel_map(pts_sf, "vis_km", "(d) Visibility (km)", "Visibility (km)") + 
#           scale_fill_viridis_c(limits = c(4, 8), option = "viridis")


#----------------------------------------------------------------------------

# ------------------ Plot composite: RH, Wind, Visibility (bubble maps) ------------------

# ── Packages ─────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(lubridate); library(purrr)
  library(sf); library(ggplot2); library(ggrepel)
  library(maps); library(patchwork); library(scales)
})

# ── Inputs ───────────────────────────────────────────────────────
# Required in memory: dust_events_cv5 with at least
# station, local_time (POSIXct), wxcodes, metar, tmpc,
# RH_CALC, wind_speed_mps, VISIBILITY_KM
CV_SHP <- "central_valley_alluvial_boundary_usgs.shp"   # <- update path if needed

# Ensure the required columns are present
required_cols <- c("station","local_time","wxcodes","metar",
                   "tmpc","RH_CALC","wind_speed_mps","VISIBILITY_KM")
missing_cols <- setdiff(required_cols, names(dust_events_cv5))
if (length(missing_cols)) {
  stop(sprintf("Missing required columns in dust_events_cv5: %s",
               paste(missing_cols, collapse=", ")))
}

# 15 stations (no FCH)
stations_df <- tibble::tribble(
  ~station, ~lat,      ~lon,
  "BFL",    35.43440, -119.0542,
  "BAB",    39.13609, -121.4366,
  "FAT",    36.78000, -119.7194,
  "HJO",    36.31139, -119.6232,
  "MAE",    36.98486, -120.1107,
  "MYV",    39.10203, -121.5688,
  "MCE",    37.28603, -120.5179,
  "OVE",    39.49000, -121.6200,
  "RBL",    40.15190, -122.2536,
  "RDD",    40.50900, -122.2934,
  "SAC",    38.50690, -121.4950,
  "SMF",    38.69542, -121.5908,
  "SCK",    37.89417, -121.2383,
  "VIS",    36.31867, -119.3929,
  "PTV",    36.02732, -119.0629
)

# ── Helpers ──────────────────────────────────────────────────────
PWC_CODES <- c("DU","BLDU","VCBLDU","DS","SS","HZ")
rx_codes   <- regex(paste0("\\b(", paste(PWC_CODES, collapse="|"), ")\\b"),
                    ignore_case = TRUE)

norm_txt  <- function(x) x %>% str_to_upper() %>% str_replace_all("VC\\s*BLDU", "VCBLDU")
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# ── 1) Identify dust days (any PWC code present that date) ───────
dust_days <- dust_events_cv5 %>%
  transmute(
    station,
    date = as.Date(local_time),
    txt  = norm_txt(paste(coalesce(as.character(wxcodes), ""),
                          coalesce(as.character(metar),   "")))
  ) %>%
  mutate(any_pwc = str_detect(txt, rx_codes)) %>%
  group_by(station, date) %>%
  summarise(dust_day = any(any_pwc, na.rm = TRUE), .groups = "drop") %>%
  filter(dust_day) %>%
  semi_join(stations_df, by = "station")

# ── 2) Daily means on dust days (per station–date) ───────────────
# *** Uses ONLY RH_CALC and VISIBILITY_KM as requested ***
daily_means <- dust_events_cv5 %>%
  mutate(date = as.Date(local_time)) %>%
  inner_join(dust_days, by = c("station","date")) %>%
  mutate(
    temp_c = suppressWarnings(as.numeric(tmpc)),
    rh_pct = suppressWarnings(as.numeric(RH_calc)),
    wind_m = suppressWarnings(as.numeric(wind_speed_mps)),
    vis_km = suppressWarnings(as.numeric(visibility_km))
  ) %>%
  group_by(station, date) %>%
  summarise(
    temp_c_day = safe_mean(temp_c),
    rh_day     = safe_mean(rh_pct),
    wind_day   = safe_mean(wind_m),
    vis_day    = safe_mean(vis_km),
    .groups = "drop"
  )

# ── 3) Composite (mean across dust days) per station ─────────────
comp_stats <- daily_means %>%
  group_by(station) %>%
  summarise(
    n_dust_days = n(),
    temp_c      = safe_mean(temp_c_day),
    rh_pct      = safe_mean(rh_day),
    wind_mps    = safe_mean(wind_day),
    vis_km      = safe_mean(vis_day),
    .groups = "drop"
  ) %>%
  right_join(stations_df, by = "station")   # keep all 15 stations (NA if no dust days)

# ── 4) Geometry: CV boundary & counties for context ──────────────
cv_boundary <- st_read(CV_SHP, quiet = TRUE) %>%
  st_zm(drop = TRUE) %>%
  st_transform(4326)

central_valley_counties <- c(
  "shasta","tehama","glenn","butte","colusa",
  "sutter","yuba","yolo","sacramento","san joaquin",
  "stanislaus","merced","madera","fresno","kings",
  "tulare","kern"
)
ca_counties <- map_data("county", "california") %>%
  filter(subregion %in% central_valley_counties)

pts_sf <- st_as_sf(comp_stats, coords = c("lon","lat"), crs = 4326)

# ── 5) Bubble-map panel (size only, no colorbar) ─────────────────
bubble_map <- function(sf_pts, value_col, title_str, legend_title,
                       size_range = c(3, 12), clamp_q = c(0.05, 0.95),
                       point_color = "#3C6E71") {
  
  vals <- sf_pts[[value_col]]
  vfin <- vals[is.finite(vals)]
  if (!length(vfin)) {
    vmin <- 0; vmax <- 1; brks <- c(0, 0.5, 1)
  } else {
    qq   <- quantile(vfin, probs = clamp_q, na.rm = TRUE, names = FALSE)
    vmin <- qq[1]; vmax <- qq[2]
    if (!is.finite(vmin) || !is.finite(vmax) || vmin == vmax) {
      vmin <- min(vfin, na.rm = TRUE); vmax <- max(vfin, na.rm = TRUE)
    }
    brks <- pretty(c(vmin, vmax), n = 4)
  }
  
  show_pts  <- is.finite(vals) & !is.na(vals)
  sf_pos    <- sf_pts[ show_pts, ]
  sf_hollow <- sf_pts[!show_pts, ]
  
  lab_df <- cbind(st_drop_geometry(sf_pts),
                  as.data.frame(st_coordinates(sf_pts)))
  names(lab_df)[(ncol(lab_df)-1):ncol(lab_df)] <- c("x","y")
  
  ggplot() +
    geom_sf(data = cv_boundary, fill = "#f0f0f0", color = "black", linewidth = 0.5) +
    geom_polygon(data = ca_counties, aes(long, lat, group = group),
                 fill = NA, color = "white", linewidth = 0.3) +
    geom_sf(data = sf_hollow, shape = 21, fill = "white", color = "grey40",
            size = size_range[1]) +
    geom_sf(
      data = sf_pos,
      aes(size = pmin(pmax(.data[[value_col]], vmin), vmax)),
      shape = 21, fill = point_color, color = "black", alpha = 0.9
    ) +
    scale_size_area(
      name = legend_title, max_size = size_range[2],
      breaks = brks, labels = scales::number_format(accuracy = 0.1)
    ) +
    ggrepel::geom_text_repel(
      data = lab_df, aes(x = x, y = y, label = station),
      size = 3, min.segment.length = 0.05, max.overlaps = 100
    ) +
    coord_sf(
      xlim = st_bbox(cv_boundary)[c("xmin","xmax")],
      ylim = st_bbox(cv_boundary)[c("ymin","ymax")],
      expand = FALSE
    ) +
    labs(title = title_str) +
    theme_void(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# ── 6) Build the three panels (size encodes value) ───────────────
p_rh <- bubble_map(
  pts_sf, "rh_pct", "(b) RH (%)", "RH (%)",
  size_range = c(3, 12), point_color = "#7B1FA2"  # purple
)

p_wind <- bubble_map(
  pts_sf, "wind_mps", "(c) Wind speed (m/s)", "Wind (m/s)",
  size_range = c(3, 12), point_color = "#D1495B"  # red
)

p_vis <- bubble_map(
  pts_sf, "vis_km", "(d) Visibility (km)", "Visibility (km)",
  size_range = c(3, 12), point_color = "#2A9D8F"  # teal
)

final_plot <- (p_rh | p_wind | p_vis)
print(final_plot)

# Optionally save:
ggsave("cv_dust_composite_bubbles.png", final_plot, width = 14, height = 6, dpi = 300)


#-----------------------------------------

# ------------------ Plot composite: RH, Wind, Visibility (discrete bubble sizes) ------------------

# ── Packages ─────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(lubridate); library(purrr)
  library(sf); library(ggplot2); library(ggrepel)
  library(maps); library(patchwork); library(scales)
})

# ── Inputs ───────────────────────────────────────────────────────
# Required in memory: dust_events_cv5 with at least
# station, local_time (POSIXct), wxcodes, metar, tmpc,
# RH_CALC, wind_speed_mps, VISIBILITY_KM
CV_SHP <- "central_valley_alluvial_boundary_usgs.shp"   # <- update path if needed

# Ensure the required columns are present
required_cols <- c("station","local_time","wxcodes","metar",
                   "tmpc","RH_CALC","wind_speed_mps","VISIBILITY_KM")
missing_cols <- setdiff(required_cols, names(dust_events_cv5))
if (length(missing_cols)) {
  stop(sprintf("Missing required columns in dust_events_cv5: %s",
               paste(missing_cols, collapse=", ")))
}

# 15 stations (no FCH)
stations_df <- tibble::tribble(
  ~station, ~lat,      ~lon,
  "BFL",    35.43440, -119.0542,
  "BAB",    39.13609, -121.4366,
  "FAT",    36.78000, -119.7194,
  "HJO",    36.31139, -119.6232,
  "MAE",    36.98486, -120.1107,
  "MYV",    39.10203, -121.5688,
  "MCE",    37.28603, -120.5179,
  "OVE",    39.49000, -121.6200,
  "RBL",    40.15190, -122.2536,
  "RDD",    40.50900, -122.2934,
  "SAC",    38.50690, -121.4950,
  "SMF",    38.69542, -121.5908,
  "SCK",    37.89417, -121.2383,
  "VIS",    36.31867, -119.3929,
  "PTV",    36.02732, -119.0629
)

# ── Helpers ──────────────────────────────────────────────────────
PWC_CODES <- c("DU","BLDU","VCBLDU","DS","SS","HZ")
rx_codes   <- regex(paste0("\\b(", paste(PWC_CODES, collapse="|"), ")\\b"),
                    ignore_case = TRUE)

norm_txt  <- function(x) x %>% str_to_upper() %>% str_replace_all("VC\\s*BLDU", "VCBLDU")
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# ── 1) Identify dust days (any PWC code present that date) ───────
dust_days <- dust_events_cv5 %>%
  transmute(
    station,
    date = as.Date(local_time),
    txt  = norm_txt(paste(coalesce(as.character(wxcodes), ""),
                          coalesce(as.character(metar),   "")))
  ) %>%
  mutate(any_pwc = str_detect(txt, rx_codes)) %>%
  group_by(station, date) %>%
  summarise(dust_day = any(any_pwc, na.rm = TRUE), .groups = "drop") %>%
  filter(dust_day) %>%
  semi_join(stations_df, by = "station")

# ── 2) Daily means on dust days (per station–date) ───────────────
# *** Uses ONLY RH_CALC and VISIBILITY_KM as requested ***
daily_means <- dust_events_cv5 %>%
  mutate(date = as.Date(local_time)) %>%
  inner_join(dust_days, by = c("station","date")) %>%
  mutate(
    temp_c = suppressWarnings(as.numeric(tmpc)),
    rh_pct = suppressWarnings(as.numeric(RH_CALC)),
    wind_m = suppressWarnings(as.numeric(wind_speed_mps)),
    vis_km = suppressWarnings(as.numeric(VISIBILITY_KM))
  ) %>%
  group_by(station, date) %>%
  summarise(
    temp_c_day = safe_mean(temp_c),
    rh_day     = safe_mean(rh_pct),
    wind_day   = safe_mean(wind_m),
    vis_day    = safe_mean(vis_km),
    .groups = "drop"
  )

# ── 3) Composite (mean across dust days) per station ─────────────
comp_stats <- daily_means %>%
  group_by(station) %>%
  summarise(
    n_dust_days = n(),
    temp_c      = safe_mean(temp_c_day),
    rh_pct      = safe_mean(rh_day),
    wind_mps    = safe_mean(wind_day),
    vis_km      = safe_mean(vis_day),
    .groups = "drop"
  ) %>%
  right_join(stations_df, by = "station")   # keep all 15 stations (NA if no dust days)

# ── 4) Geometry: CV boundary & counties for context ──────────────
cv_boundary <- st_read(CV_SHP, quiet = TRUE) %>%
  st_zm(drop = TRUE) %>%
  st_transform(4326)

central_valley_counties <- c(
  "shasta","tehama","glenn","butte","colusa",
  "sutter","yuba","yolo","sacramento","san joaquin",
  "stanislaus","merced","madera","fresno","kings",
  "tulare","kern"
)
ca_counties <- map_data("county", "california") %>%
  filter(subregion %in% central_valley_counties)

pts_sf <- st_as_sf(comp_stats, coords = c("lon","lat"), crs = 4326)

# ── 5) Discrete bubble-map panel (4 bins: S / M / L / XL) ────────
# breaks: numeric vector of length 5 (min < ... < max)
# labels_fmt: character vector of 4 labels to show in legend (e.g., "26.7–33.0")
bubble_map_discrete <- function(sf_pts, value_col, title_str, legend_title,
                                breaks, labels_fmt,
                                size_values = c(4, 7, 10, 13),
                                point_color = "#3C6E71") {
  
  stopifnot(length(breaks) == 5, length(labels_fmt) == 4, length(size_values) == 4)
  
  vals <- sf_pts[[value_col]]
  show_pts  <- is.finite(vals) & !is.na(vals)
  sf_pos    <- sf_pts[ show_pts, ]
  sf_hollow <- sf_pts[!show_pts, ]
  
  # Bin into 4 ordered categories
  words <- c("Small","Medium","Large","Very large")
  bin_cat <- cut(sf_pos[[value_col]],
                 breaks = breaks, include.lowest = TRUE, right = TRUE,
                 labels = words, ordered_result = TRUE)
  sf_pos$bin_cat <- bin_cat
  
  # Map sizes to categories
  names(size_values) <- words
  
  # For labels (stable)
  lab_df <- cbind(st_drop_geometry(sf_pts),
                  as.data.frame(st_coordinates(sf_pts)))
  names(lab_df)[(ncol(lab_df)-1):ncol(lab_df)] <- c("x","y")
  
  ggplot() +
    geom_sf(data = cv_boundary, fill = "#f0f0f0", color = "black", linewidth = 0.5) +
    geom_polygon(data = ca_counties, aes(long, lat, group = group),
                 fill = NA, color = "white", linewidth = 0.3) +
    
    # stations with NA: hollow
    geom_sf(data = sf_hollow, shape = 21, fill = "white", color = "grey40", size = 4) +
    
    # stations with values: size by binned category
    geom_sf(data = sf_pos, aes(size = bin_cat),
            shape = 21, fill = point_color, color = "black", alpha = 0.9) +
    
    scale_size_manual(
      name   = legend_title,
      values = size_values,
      breaks = words,
      labels = labels_fmt,
      guide  = guide_legend(override.aes = list(shape = 21, fill = point_color))
    ) +
    
    ggrepel::geom_text_repel(
      data = lab_df, aes(x = x, y = y, label = station),
      size = 3, min.segment.length = 0.05, max.overlaps = 100
    ) +
    coord_sf(
      xlim = st_bbox(cv_boundary)[c("xmin","xmax")],
      ylim = st_bbox(cv_boundary)[c("ymin","ymax")],
      expand = FALSE
    ) +
    labs(title = title_str) +
    theme_void(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# Helper to make equal-interval breaks & "a–b" labels with chosen precision
mk_breaks_labels <- function(vmin, vmax, digits = 1) {
  b <- seq(vmin, vmax, length.out = 5)
  lab <- paste0(format(round(b[-5], digits), nsmall = digits), "–",
                format(round(b[-1], digits), nsmall = digits))
  list(breaks = b, labels = lab)
}

# ── 6) Build the three panels with YOUR fixed ranges ─────────────
# RH (%): 26.66 → 51.93 (1 decimal)
rh_bl <- mk_breaks_labels(26.66, 51.93, digits = 1)
p_rh <- bubble_map_discrete(
  pts_sf, "rh_pct", "(b) RH (%)", "RH (%)",
  breaks = rh_bl$breaks, labels_fmt = rh_bl$labels,
  size_values = c(4, 7, 10, 13), point_color = "#7B1FA2"
)

# Wind (m/s): 6.842 → 9.214 (2 decimals for clarity)
wind_bl <- mk_breaks_labels(6.842, 9.214, digits = 2)
p_wind <- bubble_map_discrete(
  pts_sf, "wind_mps", "(c) Wind speed (m/s)", "Wind (m/s)",
  breaks = wind_bl$breaks, labels_fmt = wind_bl$labels,
  size_values = c(4, 7, 10, 13), point_color = "#D1495B"
)

# Visibility (km): 6.183 → 10.294 (1 decimal)
vis_bl <- mk_breaks_labels(6.183, 10.294, digits = 1)
p_vis <- bubble_map_discrete(
  pts_sf, "vis_km", "(d) Visibility (km)", "Visibility (km)",
  breaks = vis_bl$breaks, labels_fmt = vis_bl$labels,
  size_values = c(4, 7, 10, 13), point_color = "#2A9D8F"
)

final_plot <- (p_rh | p_wind | p_vis)
print(final_plot)

# Optionally save:
ggsave("cv_dust_composite_bubbles_discrete.png", final_plot, width = 14, height = 6, dpi = 300)






#-----------------------------
# ------------------ CV dust composite: RH, Wind, Visibility (3-size bubbles, fixed bins) ------------------

# ── Packages ─────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(lubridate); library(purrr)
  library(sf); library(ggplot2); library(ggrepel)
  library(maps); library(patchwork); library(scales)
})

# ── Inputs ───────────────────────────────────────────────────────
# Required in memory: `dust_events_cv5` with at least:
# station, local_time (POSIXct), wxcodes, metar, tmpc, RH_CALC, wind_speed_mps, VISIBILITY_KM
CV_SHP <- "central_valley_alluvial_boundary_usgs.shp"   # <- update path if needed

# Ensure required columns exist
required_cols <- c("station","local_time","wxcodes","metar",
                   "tmpc","RH_CALC","wind_speed_mps","VISIBILITY_KM")
missing_cols <- setdiff(required_cols, names(dust_events_cv5))
if (length(missing_cols)) stop("Missing columns in dust_events_cv5: ", paste(missing_cols, collapse=", "))

# 15 stations (no FCH)
stations_df <- tibble::tribble(
  ~station, ~lat,      ~lon,
  "BFL",    35.43440, -119.0542,
  "BAB",    39.13609, -121.4366,
  "FAT",    36.78000, -119.7194,
  "HJO",    36.31139, -119.6232,
  "MAE",    36.98486, -120.1107,
  "MYV",    39.10203, -121.5688,
  "MCE",    37.28603, -120.5179,
  "OVE",    39.49000, -121.6200,
  "RBL",    40.15190, -122.2536,
  "RDD",    40.50900, -122.2934,
  "SAC",    38.50690, -121.4950,
  "SMF",    38.69542, -121.5908,
  "SCK",    37.89417, -121.2383,
  "VIS",    36.31867, -119.3929,
  "PTV",    36.02732, -119.0629
)

# ── Helpers ──────────────────────────────────────────────────────
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
PWC_CODES <- c("DU","BLDU","VCBLDU","DS","SS","HZ")
rx_codes   <- regex(paste0("\\b(", paste(PWC_CODES, collapse="|"), ")\\b"), ignore_case = TRUE)
norm_txt   <- function(x) x |> str_to_upper() |> str_replace_all("VC\\s*BLDU", "VCBLDU")

# ── 1) Identify dust days (any PWC code present that date) ───────
dust_days <- dust_events_cv5 %>%
  transmute(
    station,
    date = as.Date(local_time),
    txt  = norm_txt(paste(coalesce(as.character(wxcodes), ""),
                          coalesce(as.character(metar),   "")))
  ) %>%
  mutate(any_pwc = str_detect(txt, rx_codes)) %>%
  group_by(station, date) %>%
  summarise(dust_day = any(any_pwc, na.rm = TRUE), .groups = "drop") %>%
  filter(dust_day) %>%
  semi_join(stations_df, by = "station")

# ── 2) Daily means on dust days (per station–date) ───────────────
# *** ONLY RH_CALC and VISIBILITY_KM are used ***
daily_means <- dust_events_cv5 %>%
  mutate(date = as.Date(local_time)) %>%
  inner_join(dust_days, by = c("station","date")) %>%
  mutate(
    temp_c = suppressWarnings(as.numeric(tmpc)),
    rh_pct = suppressWarnings(as.numeric(RH_CALC)),
    wind_m = suppressWarnings(as.numeric(wind_speed_mps)),
    vis_km = suppressWarnings(as.numeric(VISIBILITY_KM))
  ) %>%
  group_by(station, date) %>%
  summarise(
    temp_c_day = safe_mean(temp_c),
    rh_day     = safe_mean(rh_pct),
    wind_day   = safe_mean(wind_m),
    vis_day    = safe_mean(vis_km),
    .groups = "drop"
  )

# ── 3) Composite (mean across dust days) per station ─────────────
comp_stats <- daily_means %>%
  group_by(station) %>%
  summarise(
    n_dust_days = n(),
    temp_c      = safe_mean(temp_c_day),
    rh_pct      = safe_mean(rh_day),
    wind_mps    = safe_mean(wind_day),
    vis_km      = safe_mean(vis_day),
    .groups = "drop"
  ) %>%
  right_join(stations_df, by = "station")   # keep all 15 stations (NA if no dust days)

# ── 4) Geometry: CV boundary & counties for context ──────────────
cv_boundary <- st_read(CV_SHP, quiet = TRUE) %>%
  st_zm(drop = TRUE) %>%
  st_transform(4326)

central_valley_counties <- c("shasta","tehama","glenn","butte","colusa",
                             "sutter","yuba","yolo","sacramento","san joaquin",
                             "stanislaus","merced","madera","fresno","kings",
                             "tulare","kern")
ca_counties <- map_data("county", "california") %>%
  filter(subregion %in% central_valley_counties)

pts_sf <- st_as_sf(comp_stats, coords = c("lon","lat"), crs = 4326)

# ── 5) Generic discrete bubble map (N bins) ──────────────────────
bubble_map_bins <- function(sf_pts, value_col, title_str, legend_title,
                            breaks, labels_fmt,
                            size_values = c(6, 10, 14),
                            point_color = "#3C6E71") {
  # breaks: numeric of length N+1; labels_fmt: length N; size_values: length N
  N <- length(breaks) - 1
  stopifnot(N >= 1, length(labels_fmt) == N, length(size_values) == N)
  
  vals <- sf_pts[[value_col]]
  show_pts  <- is.finite(vals) & !is.na(vals)
  sf_pos    <- sf_pts[ show_pts, ]
  sf_hollow <- sf_pts[!show_pts, ]
  
  # include max by extending top break slightly
  breaks[length(breaks)] <- breaks[length(breaks)] + 1e-9
  lvl_names <- paste0("bin", seq_len(N))
  sf_pos$bin_cat <- cut(sf_pos[[value_col]],
                        breaks = breaks, include.lowest = TRUE, right = TRUE,
                        labels = lvl_names, ordered_result = TRUE)
  
  # Labels/coords
  lab_df <- cbind(st_drop_geometry(sf_pts), as.data.frame(st_coordinates(sf_pts)))
  names(lab_df)[(ncol(lab_df)-1):ncol(lab_df)] <- c("x","y")
  
  ggplot() +
    geom_sf(data = cv_boundary, fill = "#f0f0f0", color = "black", linewidth = 0.5) +
    geom_polygon(data = ca_counties, aes(long, lat, group = group),
                 fill = NA, color = "white", linewidth = 0.3) +
    geom_sf(data = sf_hollow, shape = 21, fill = "white", color = "grey40", size = size_values[1]) +
    geom_sf(data = sf_pos, aes(size = bin_cat),
            shape = 21, fill = point_color, color = "black", alpha = 0.9) +
    scale_size_manual(
      name   = legend_title,
      values = setNames(size_values, lvl_names),
      breaks = lvl_names,
      labels = labels_fmt,
      drop   = FALSE,
      guide  = guide_legend(override.aes = list(shape = 21, fill = point_color))
    ) +
    ggrepel::geom_text_repel(
      data = lab_df, aes(x = x, y = y, label = station),
      size = 3, min.segment.length = 0.05, max.overlaps = 100
    ) +
    coord_sf(
      xlim = st_bbox(cv_boundary)[c("xmin","xmax")],
      ylim = st_bbox(cv_boundary)[c("ymin","ymax")],
      expand = FALSE
    ) +
    labs(title = title_str) +
    theme_void(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# ── 6) Three-bin settings (YOUR EXACT BOUNDS) --------------------
# RH (%): 0–30, 31–45, 46–55
rh_breaks <- c(0, 30, 45, 55)
rh_labels <- c("0–30", "31–45", "46–55")

# Wind (m/s): 6–7.5, 7.6–8.6, 8.7–10
wind_breaks <- c(6.0, 7.5, 8.6, 10.0)
wind_labels <- c("6–7.5", "7.6–8.6", "8.7–10")

# Visibility (km): 6–7.2, 7.3–8.2, 8.3–11
vis_breaks <- c(6.0, 7.2, 8.2, 11.0)     # results in [6,7.2], (7.2,8.2], (8.2,11]
vis_labels <- c("6–7.2", "7.3–8.2", "8.3–11")

# ── 7) Build the three panels -----------------------------------
p_rh <- bubble_map_bins(
  pts_sf, "rh_pct", "(c) Relative humidity", "RH (%)",
  breaks = rh_breaks, labels_fmt = rh_labels,
  size_values = c(6, 10, 14), point_color = "#7B1FA2" # purple
)

p_wind <- bubble_map_bins(
  pts_sf, "wind_mps", "(a) Wind speed", "Wind (m/s)",
  breaks = wind_breaks, labels_fmt = wind_labels,
  size_values = c(6, 10, 14), point_color = "#D1495B" # red
)

p_vis <- bubble_map_bins(
  pts_sf, "vis_km", "(b) Visibility", "Visibility (km)",
  breaks = vis_breaks, labels_fmt = vis_labels,
  size_values = c(6, 10, 14), point_color = "#2A9D8F" # teal
)

final_plot <- (p_rh | p_wind | p_vis)
final_plot <- (p_wind | p_vis | p_rh)
print(final_plot)

# Optionally save:
ggsave("cv_dust_composite_bubbles_3sizes_fixed.png", final_plot, width = 14, height = 6, dpi = 300)











#----------------------------- dust algorithm lighting plot 

# Load necessary libraries
library(ggplot2)
library(maps)
library(dplyr)

# Define Central Valley counties
central_valley_counties <- c("shasta", "tehama", "glenn", "butte", "colusa",
                             "sutter", "yuba", "yolo", "sacramento", "san joaquin",
                             "stanislaus", "merced", "madera", "fresno", "kings",
                             "tulare", "kern")

# Get map data for California counties and filter for Central Valley counties
california_counties <- map_data("county", region = "california")
central_valley_map <- california_counties %>% filter(subregion %in% central_valley_counties)

# Get California state map data
california_map <- map_data("state", region = "california")

library(dplyr)
library(lubridate)
library(riem)
library(purrr)


# ── 3) Build the station coordinate table ────────────────────────
df_coords <- tribble(
  ~station, ~lat,      ~lon,
  "BFL",    35.43440, -119.0542,
  "BAB",    39.13609, -121.4366,
  "FAT",    36.78000, -119.7194,
  "HJO",    36.31139, -119.6232,
  "MAE",    36.98486, -120.1107,
  "MYV",    39.10203, -121.5688,
  "MCE",    37.28603, -120.5179,
  "OVE",    39.49000, -121.6200,
  "PTV",    36.02732, -119.0629,
  "RBL",    40.15190, -122.2536,
  "RDD",    40.50900, -122.2934,
  "SAC",    38.50690, -121.4950,
  "SMF",    38.69542, -121.5908,
  "SCK",    37.89417, -121.2383,
  "VIS",    36.31867, -119.3929
)


stations_sf <- st_as_sf(df_coords, coords = c("lon","lat"), crs = 4326)

#── 3) Read & reproject USGS Central Valley alluvial boundary ──────
cv_boundary <- st_read("central_valley_alluvial_boundary_usgs.shp", quiet = TRUE) %>%
  st_zm(drop = TRUE) %>%
  st_transform(4326)

#── 4) Draw Central Valley county lines via maps package ───────────
central_valley_counties <- c("shasta","tehama","glenn","butte","colusa",
                             "sutter","yuba","yolo","sacramento","san joaquin",
                             "stanislaus","merced","madera","fresno","kings",
                             "tulare","kern")
ca_counties <- map_data("county", "california") %>%
  filter(subregion %in% central_valley_counties)

# 1) fetch each station separately and row-bind
df_raw <- map_dfr(df_coords$station, function(st) {
  riem_measures(
    station    = st,
    date_start = "2021-10-11",
    date_end   = "2021-10-11"
  )
})

# 2) keep only hours 18–21 UTC, convert / round, join coords, and format
df_table <- df_raw %>%
  filter(hour(valid) %in% 18:21) %>%
  transmute(
    station,
    Category          = paste0(sprintf("%02d", hour(valid)), ":00UTC"),
    `Visibility` = round(vsby * 1.60934, 2),
    `Wind-speed`= round(sknt * 0.514444, 2),
    `Temperature`= round((tmpf - 32) * 5/9, 1)
  ) %>%
  left_join(df_coords, by = "station") %>%
  arrange(station, Category)

print(df_table)

str(df_table)

df_table$Category <- factor(df_table$Category, levels = c("18:00UTC", "19:00UTC", "20:00UTC","21:00UTC"))

df_table$Category <- factor(
  df_table$Category,
  levels = c("18:00UTC", "19:00UTC", "20:00UTC", "21:00UTC"),
  labels = c("11:00 LST", "12:00 LST", "13:00 LST", "14:00 LST")
)


# Define Central Valley counties
central_valley_counties <- c("shasta", "tehama", "glenn", "butte", "colusa",
                             "sutter", "yuba", "yolo", "sacramento", "san joaquin",
                             "stanislaus", "merced", "madera", "fresno", "kings",
                             "tulare", "kern")

# Get map data for California counties and filter for Central Valley counties
california_counties <- map_data("county", region = "california")
central_valley_map <- california_counties %>% filter(subregion %in% central_valley_counties)

# Get California state map data
california_map <- map_data("state", region = "california")

png("tdust_plot.png", width = 1150, height = 350)

# Create a three-panel grid
png("VIS_plot.png", width = 1150, height = 550)
ggplot() +
  geom_sf(data = cv_boundary,
          fill  = "#f0f0f0",
          color = "black",
          size  = 7) +
  # b) county subdivisions
  geom_polygon(data = ca_counties,
               aes(x = long, y = lat, group = group),
               fill  = NA,
               color = "white",
               size  = 0.5) +
  geom_polygon(data = california_map, aes(x = long, y = lat, group = group),
               fill = "white", color = "black", alpha = 0.4) +
  # Plot dust data points with color representing Dust count and size based on Dust
  geom_point(data = df_table, aes(x = lon, y = lat, color = Visibility, size = 0.8), alpha = 0.8) +
  # — overlay an “x” on any point with visibility < 10 km —
  geom_point(data = subset(df_table, Visibility < 10 & `Wind-speed` >= 6),
             aes(x = lon, y = lat),
             shape = 4,        # ×
             size = 3,         # adjust to taste
             stroke = 1.2,     # line thickness
             color = "black") +
  # Set a multi-step color gradient for Dust count
  scale_color_gradientn(colors = c("yellow", "orange", "red", "darkred"),
                        limits = c(0, 25), breaks = seq(0, 25, by = 5),
                        name = "Visibility (Km)") +
  # Add title
  labs(title = " ") +
  # Adjust theme to remove gridlines, axis ticks, and axis labels
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",               # Position legend to the right (vertical orientation)
    panel.grid.major = element_blank(),       # Remove major gridlines
    panel.grid.minor = element_blank(),       # Remove minor gridlines
    axis.title = element_blank(),             # Remove axis titles
    axis.text = element_blank(),              # Remove axis text (tick labels)
    axis.ticks = element_blank()              # Remove axis ticks
  ) +
  guides(
    color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barwidth = unit(0.4, "cm"),
                           barheight = unit(8, "cm")),
    size = "none"                             # Remove the size legend
  ) +
  # Create a three-panel grid
  facet_wrap(~ Category, ncol = 4)
dev.off()





ggplot() +
  geom_sf(data = cv_boundary,
          fill  = "#f0f0f0",
          color = "black",
          size  = 7) +
  # b) county subdivisions
  geom_polygon(data = ca_counties,
               aes(x = long, y = lat, group = group),
               fill  = NA,
               color = "white",
               size  = 0.5) +
  geom_polygon(data = california_map, aes(x = long, y = lat, group = group),
               fill = "white", color = "black", alpha = 0.4) +
  # Plot dust data points with color representing Dust count and size based on Dust
  geom_point(data = df_table, aes(x = lon, y = lat, color = `Wind-speed`, size = 0.8), alpha = 0.8) +
  # — overlay an “x” on any point with visibility < 10 km —
  geom_point(data = subset(df_table, Visibility < 10 & `Wind-speed` >= 6),
             aes(x = lon, y = lat),
             shape = 4,        # ×
             size = 3,         # adjust to taste
             stroke = 1.2,     # line thickness
             color = "black") +
  # Set a multi-step color gradient for Dust count
  scale_color_gradientn(colors = c("yellow", "orange", "red", "darkred"),
                        limits = c(0, 25), breaks = seq(0, 25, by = 5),
                        name = "Wind-Speed (m/s)") +
  # Add title
  labs(title = " ") +
  # Adjust theme to remove gridlines, axis ticks, and axis labels
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",               # Position legend to the right (vertical orientation)
    panel.grid.major = element_blank(),       # Remove major gridlines
    panel.grid.minor = element_blank(),       # Remove minor gridlines
    axis.title = element_blank(),             # Remove axis titles
    axis.text = element_blank(),              # Remove axis text (tick labels)
    axis.ticks = element_blank()              # Remove axis ticks
  ) +
  
  guides(
    color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barwidth = unit(0.4, "cm"),
                           barheight = unit(8, "cm")),
    size = "none"                             # Remove the size legend
  ) +
  # Create a three-panel grid
  facet_wrap(~ Category, ncol = 4)

#-------------------------------------END----------------------------------------------

#---------------------------------PLOT FOR 2024

#----------------------------- dust algorithm lighting plot 

# Load necessary libraries
library(ggplot2)
library(maps)
library(dplyr)

# Define Central Valley counties
central_valley_counties <- c("shasta", "tehama", "glenn", "butte", "colusa",
                             "sutter", "yuba", "yolo", "sacramento", "san joaquin",
                             "stanislaus", "merced", "madera", "fresno", "kings",
                             "tulare", "kern")

# Get map data for California counties and filter for Central Valley counties
california_counties <- map_data("county", region = "california")
central_valley_map <- california_counties %>% filter(subregion %in% central_valley_counties)

# Get California state map data
california_map <- map_data("state", region = "california")

library(dplyr)
library(lubridate)
library(riem)
library(purrr)


# ── 3) Build the station coordinate table ────────────────────────
df_coords <- tribble(
  ~station, ~lat,      ~lon,
  "BFL",    35.43440, -119.0542,
  "BAB",    39.13609, -121.4366,
  "FAT",    36.78000, -119.7194,
  "HJO",    36.31139, -119.6232,
  "MAE",    36.98486, -120.1107,
  "MYV",    39.10203, -121.5688,
  "MCE",    37.28603, -120.5179,
  "OVE",    39.49000, -121.6200,
  "PTV",    36.02732, -119.0629,
  "RBL",    40.15190, -122.2536,
  "RDD",    40.50900, -122.2934,
  "SAC",    38.50690, -121.4950,
  "SMF",    38.69542, -121.5908,
  "SCK",    37.89417, -121.2383,
  "VIS",    36.31867, -119.3929
)


stations_sf <- st_as_sf(df_coords, coords = c("lon","lat"), crs = 4326)

#── 3) Read & reproject USGS Central Valley alluvial boundary ──────
cv_boundary <- st_read("central_valley_alluvial_boundary_usgs.shp", quiet = TRUE) %>%
  st_zm(drop = TRUE) %>%
  st_transform(4326)

#── 4) Draw Central Valley county lines via maps package ───────────
central_valley_counties <- c("shasta","tehama","glenn","butte","colusa",
                             "sutter","yuba","yolo","sacramento","san joaquin",
                             "stanislaus","merced","madera","fresno","kings",
                             "tulare","kern")
ca_counties <- map_data("county", "california") %>%
  filter(subregion %in% central_valley_counties)

# 1) fetch each station separately and row-bind
df_raw <- map_dfr(df_coords$station, function(st) {
  riem_measures(
    station    = st,
    date_start = "2024-11-11",
    date_end   = "2024-11-11"
  )
})

# 2) keep only hours 18–21 UTC, convert / round, join coords, and format
df_table <- df_raw %>%
  filter(hour(valid) %in% 20:23) %>%
  transmute(
    station,
    Category          = paste0(sprintf("%02d", hour(valid)), ":00UTC"),
    `Visibility` = round(vsby * 1.60934, 2),
    `Wind-speed`= round(sknt * 0.514444, 2),
    `Temperature`= round((tmpf - 32) * 5/9, 1)
  ) %>%
  left_join(df_coords, by = "station") %>%
  arrange(station, Category)

print(df_table)

str(df_table)

df_table$Category <- factor(df_table$Category, levels = c("20:00UTC", "21:00UTC", "22:00UTC","23:00UTC"))

df_table$Category <- factor(
  df_table$Category,
  levels = c("20:00UTC", "21:00UTC", "22:00UTC", "23:00UTC"),
  labels = c("13:00 LST", "14:00 LST", "15:00 LST", "16:00 LST")
)


# Define Central Valley counties
central_valley_counties <- c("shasta", "tehama", "glenn", "butte", "colusa",
                             "sutter", "yuba", "yolo", "sacramento", "san joaquin",
                             "stanislaus", "merced", "madera", "fresno", "kings",
                             "tulare", "kern")

# Get map data for California counties and filter for Central Valley counties
california_counties <- map_data("county", region = "california")
central_valley_map <- california_counties %>% filter(subregion %in% central_valley_counties)

# Get California state map data
california_map <- map_data("state", region = "california")

png("tdust_plot.png", width = 1150, height = 350)

# Create a three-panel grid
png("VIS_plot.png", width = 1150, height = 550)
ggplot() +
  geom_sf(data = cv_boundary,
          fill  = "#f0f0f0",
          color = "black",
          size  = 7) +
  # b) county subdivisions
  geom_polygon(data = ca_counties,
               aes(x = long, y = lat, group = group),
               fill  = NA,
               color = "white",
               size  = 0.5) +
  geom_polygon(data = california_map, aes(x = long, y = lat, group = group),
               fill = "white", color = "black", alpha = 0.4) +
  # Plot dust data points with color representing Dust count and size based on Dust
  geom_point(data = df_table, aes(x = lon, y = lat, color = Visibility, size = 0.8), alpha = 0.8) +
  # — overlay an “x” on any point with visibility < 10 km —
  geom_point(data = subset(df_table, Visibility < 10 & `Wind-speed` >= 6),
             aes(x = lon, y = lat),
             shape = 4,        # ×
             size = 3,         # adjust to taste
             stroke = 1.2,     # line thickness
             color = "black") +
  # Set a multi-step color gradient for Dust count
  scale_color_gradientn(colors = c("yellow", "orange", "red", "darkred"),
                        limits = c(0, 25), breaks = seq(0, 25, by = 5),
                        name = "Visibility (Km)") +
  # Add title
  labs(title = " ") +
  # Adjust theme to remove gridlines, axis ticks, and axis labels
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",               # Position legend to the right (vertical orientation)
    panel.grid.major = element_blank(),       # Remove major gridlines
    panel.grid.minor = element_blank(),       # Remove minor gridlines
    axis.title = element_blank(),             # Remove axis titles
    axis.text = element_blank(),              # Remove axis text (tick labels)
    axis.ticks = element_blank()              # Remove axis ticks
  ) +
  guides(
    color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barwidth = unit(0.4, "cm"),
                           barheight = unit(8, "cm")),
    size = "none"                             # Remove the size legend
  ) +
  # Create a three-panel grid
  facet_wrap(~ Category, ncol = 4)
dev.off()





ggplot() +
  geom_sf(data = cv_boundary,
          fill  = "#f0f0f0",
          color = "black",
          size  = 7) +
  # b) county subdivisions
  geom_polygon(data = ca_counties,
               aes(x = long, y = lat, group = group),
               fill  = NA,
               color = "white",
               size  = 0.5) +
  geom_polygon(data = california_map, aes(x = long, y = lat, group = group),
               fill = "white", color = "black", alpha = 0.4) +
  # Plot dust data points with color representing Dust count and size based on Dust
  geom_point(data = df_table, aes(x = lon, y = lat, color = `Wind-speed`, size = 0.8), alpha = 0.8) +
  # — overlay an “x” on any point with visibility < 10 km —
  geom_point(data = subset(df_table, Visibility < 10 & `Wind-speed` >= 6),
             aes(x = lon, y = lat),
             shape = 4,        # ×
             size = 3,         # adjust to taste
             stroke = 1.2,     # line thickness
             color = "black") +
  # Set a multi-step color gradient for Dust count
  scale_color_gradientn(colors = c("yellow", "orange", "red", "darkred"),
                        limits = c(0, 25), breaks = seq(0, 25, by = 5),
                        name = "Wind-Speed (m/s)") +
  # Add title
  labs(title = " ") +
  # Adjust theme to remove gridlines, axis ticks, and axis labels
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",               # Position legend to the right (vertical orientation)
    panel.grid.major = element_blank(),       # Remove major gridlines
    panel.grid.minor = element_blank(),       # Remove minor gridlines
    axis.title = element_blank(),             # Remove axis titles
    axis.text = element_blank(),              # Remove axis text (tick labels)
    axis.ticks = element_blank()              # Remove axis ticks
  ) +
  
  guides(
    color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barwidth = unit(0.4, "cm"),
                           barheight = unit(8, "cm")),
    size = "none"                             # Remove the size legend
  ) +
  # Create a three-panel grid
  facet_wrap(~ Category, ncol = 4)

#---------------------------------------------------------------




#-------------------------climatology----------
# -------------------- Packages --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(forcats)
library(patchwork)

# -------------------- Config ----------------------
stations_15 <- c("BFL","BAB","FAT","HJO","MAE","MYV","MCE","OVE",
                 "RBL","RDD","SAC","SMF","SCK","VIS","PTV")

# Helper to add time fields and coerce numerics
prep_time <- function(df) {
  df %>%
    mutate(
      year  = year(local_time),
      month = factor(month.abb[month(local_time)],
                     levels = month.abb, ordered = TRUE),
      wind  = suppressWarnings(as.numeric(wind_speed_mps)),
      vis   = suppressWarnings(as.numeric(visibility_km)),
      rh    = suppressWarnings(as.numeric(RH_calc))
    )
}

# -------------------- Data prep --------------------
# Dust dataset (already dust-only rows)
dust <- dust_events_cv5 %>%
  filter(station %in% stations_15) %>%
  prep_time()

# Year window based on dust data
yr_min <- min(dust$year, na.rm = TRUE)
yr_max <- max(dust$year, na.rm = TRUE)

# Climatology (all rows), restricted to same stations + years,
# and FILTER OUT calm (0 m/s) hours, but KEEP NA winds.
climo <- asos_data9 %>%
  filter(station %in% stations_15,
         between(year(local_time), yr_min, yr_max)) %>%
  prep_time() %>%
  filter(is.na(wind) | wind > 0)

# -------------------- Summaries --------------------
yearly_summary <- function(df, label) {
  df %>%
    group_by(year) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-year, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

monthly_summary <- function(df, label) {
  df %>%
    group_by(month) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-month, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

yearly_all  <- bind_rows(
  yearly_summary(climo, "Climatology"),
  yearly_summary(dust,  "Dust events")
)

monthly_all <- bind_rows(
  monthly_summary(climo, "Climatology"),
  monthly_summary(dust,  "Dust events")
)

# -------------------- Plotting ---------------------
pal <- c(`Climatology` = "grey80", `Dust events` = "#2C7FB8")

base_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  )

# Yearly (top row)
p_year <- ggplot(yearly_all, aes(x = year, y = value, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pal) +
  labs(y = "Yearly mean") +
  base_theme +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))

# Monthly (bottom row)
p_mon <- ggplot(monthly_all, aes(x = month, y = value, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pal) +
  labs(y = "Monthly mean", x = "Month") +
  base_theme

final_plot <- p_year / p_mon 

final_plot


pm10=read.csv('pm_data.csv')
unique(pm10$Station)


#---------------------------------------------------------------
# -------------------- Packages --------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(forcats)
library(patchwork)
library(ggpattern)   # install.packages("ggpattern") if needed

# --------------- Keep your data prep as you had ---------------
# (Same as before: build 'dust', 'climo', 'yearly_all', 'monthly_all')
# NOTE: climatology filters calm hours (wind == 0) but keeps NA

stations_15 <- c("BFL","BAB","FAT","HJO","MAE","MYV","MCE","OVE",
                 "RBL","RDD","SAC","SMF","SCK","VIS","PTV")

prep_time <- function(df) {
  df %>%
    mutate(
      year  = year(local_time),
      month = factor(month.abb[month(local_time)],
                     levels = month.abb, ordered = TRUE),
      wind  = suppressWarnings(as.numeric(wind_speed_mps)),
      vis   = suppressWarnings(as.numeric(visibility_km)),
      rh    = suppressWarnings(as.numeric(RH_calc))
    )
}

dust <- dust_events_cv5 %>%
  filter(station %in% stations_15) %>%
  prep_time()

yr_min <- min(dust$year, na.rm = TRUE)
yr_max <- max(dust$year, na.rm = TRUE)

climo <- asos_data9 %>%
  filter(station %in% stations_15,
         between(year(local_time), yr_min, yr_max)) %>%
  prep_time() %>%
  filter(is.na(wind) | wind > 0)   # <<< remove calm hours only

yearly_summary <- function(df, label) {
  df %>%
    group_by(year) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-year, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

monthly_summary <- function(df, label) {
  df %>%
    group_by(month) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-month, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

yearly_all  <- bind_rows(
  yearly_summary(climo, "Climatology"),
  yearly_summary(dust,  "Dust events")
)

monthly_all <- bind_rows(
  monthly_summary(climo, "Climatology"),
  monthly_summary(dust,  "Dust events")
)

# -------------------- Reorder facets --------------------
metric_levels <- c("Wind speed (m/s)", "Visibility (km)", "Relative humidity (%)")
yearly_all  <- yearly_all  %>% mutate(metric = fct_relevel(metric, metric_levels))
monthly_all <- monthly_all %>% mutate(metric = fct_relevel(metric, metric_levels))

# -------------------- ggpattern styling --------------------
fill_cols   <- c(`Climatology` = "grey80", `Dust events` = "#2C7FB8")
stripe_cols <- c(`Climatology` = "grey80", `Dust events` = "#1B5E9A")

base_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  )

bar_plot <- function(d, x_var, ylab, xlab = NULL) {
  ggplot(
    d,
    aes(
      x = {{x_var}},
      y = value,
      fill = dataset,
      pattern = dataset,
      pattern_fill = dataset
    )
  ) +
    ggpattern::geom_col_pattern(
      position = position_dodge(width = 0.75),
      width = 0.65,
      colour = "grey20",
      pattern_colour = "grey20",
      pattern_density = 0.35,
      pattern_spacing = 0.02,
      pattern_size = 0.25,
      pattern_angle = ifelse(d$dataset == "Dust events", 45, 0)
    ) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = fill_cols) +
    scale_pattern_manual(values = c(`Climatology` = "none", `Dust events` = "stripe")) +
    scale_pattern_fill_manual(values = stripe_cols) +
    labs(y = ylab, x = xlab) +
    base_theme
}

# -------------------- Build final figure --------------------
p_year <- bar_plot(yearly_all,  year,   ylab = "Yearly mean") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9))

p_mon  <- bar_plot(monthly_all, month,  ylab = "Monthly mean", xlab = "Month")

final_plot <- p_year / p_mon 

final_plot
ggsave("monthly_pattern.png", final_plot, width = 18, height = 8, dpi = 300)





#--------------------------PM10------------------------------------

#---------------------------------------------------------------
# Packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(forcats)
  library(patchwork)
  library(ggpattern)   # install.packages("ggpattern") if needed
  library(readr)
})

# -------------------- Station sets --------------------
stations_15 <- c("BFL","BAB","FAT","HJO","MAE","MYV","MCE","OVE",
                 "RBL","RDD","SAC","SMF","SCK","VIS","PTV")  # ASOS vars
pm_targets  <- c("HJO","BAB","BFL","FAT","MAE","MCE","MYV","SAC","SCK","VIS")  # PM network

# -------------------- Helpers --------------------
prep_time <- function(df) {
  df %>%
    mutate(
      year  = lubridate::year(local_time),
      month = factor(month.abb[lubridate::month(local_time)],
                     levels = month.abb, ordered = TRUE),
      wind  = suppressWarnings(as.numeric(wind_speed_mps)),
      vis   = suppressWarnings(as.numeric(visibility_km)),
      rh    = suppressWarnings(as.numeric(RH_calc))
    )
}

yearly_summary <- function(df, label) {
  df %>%
    group_by(year) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-year, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

monthly_summary <- function(df, label) {
  df %>%
    group_by(month) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-month, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

mk_month_ord <- function(x) factor(as.character(x), levels = month.abb, ordered = TRUE)

# -------------------- Build ASOS sets (your original logic) --------------------
# Objects 'dust_events_cv5' and 'asos_data9' must be available in your session.
dust <- dust_events_cv5 %>%
  filter(station %in% stations_15) %>%
  prep_time()

yr_min <- min(dust$year, na.rm = TRUE)
yr_max <- max(dust$year, na.rm = TRUE)

climo <- asos_data9 %>%
  filter(station %in% stations_15,
         between(lubridate::year(local_time), yr_min, yr_max)) %>%
  prep_time() %>%
  filter(is.na(wind) | wind > 0)   # remove calm hours only

yearly_all  <- bind_rows(
  yearly_summary(climo, "Climatology"),
  yearly_summary(dust,  "Dust events")
)

monthly_all <- bind_rows(
  monthly_summary(climo, "Climatology"),
  monthly_summary(dust,  "Dust events")
)

# -------------------- NEW: PM10 (climatology vs dusty days) --------------------
pm_file   <- "pm_data.csv"
dust_file <- "dusty-stations.csv"

# Load all-days PM
pm_all <- read_csv(pm_file, show_col_types = FALSE) %>%
  mutate(
    Station = toupper(trimws(Station)),
    dt      = parse_date_time(
      Date,
      orders = c("mdy HMS","mdy HM","mdy H",
                 "Ymd HMS","Y-m-d H:M:S","Y-m-d H:M",
                 "dmy HMS","dmy HM"),
      quiet  = TRUE
    ),
    date    = as_date(dt)
  ) %>%
  filter(Station %in% pm_targets, !is.na(PM10), !is.na(dt))

# Load dusty-day keys (LOCAL date)
dust_keys <- dust_events_cv5 %>%
  transmute(
    Station = toupper(trimws(station)),
    date    = as_date(parse_date_time(date, orders = c("Y-m-d","Y/m/d","mdy"),
                                      exact = FALSE, quiet = TRUE))
  ) %>%
  filter(Station %in% pm_targets) %>%
  distinct()

# PM on dusty days (join by Station + LOCAL date)
pm_dust <- pm_all %>% inner_join(dust_keys, by = c("Station","date"))

# ---- Yearly PM10 ----
yearly_pm_all <- pm_all %>%
  mutate(year = year(dt)) %>%
  group_by(year) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Climatology")

yearly_pm_dust <- pm_dust %>%
  mutate(year = year(dt)) %>%
  group_by(year) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Dust events")

# ---- Monthly PM10 climatology ----
monthly_pm_all <- pm_all %>%
  mutate(month = factor(month(dt), levels = 1:12, labels = month.abb, ordered = TRUE)) %>%
  group_by(month) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Climatology")

monthly_pm_dust <- pm_dust %>%
  mutate(month = factor(month(dt), levels = 1:12, labels = month.abb, ordered = TRUE)) %>%
  group_by(month) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Dust events")

# ---- Harmonize factor types (fix bind_rows error) ----
monthly_pm_all  <- monthly_pm_all  %>% mutate(month = mk_month_ord(month))
monthly_pm_dust <- monthly_pm_dust %>% mutate(month = mk_month_ord(month))
monthly_all     <- monthly_all     %>% mutate(month = mk_month_ord(month))

# ---- Append PM10 to the ASOS long tables ----
yearly_all  <- bind_rows(
  yearly_all,
  yearly_pm_all  %>% select(year, metric, value, dataset),
  yearly_pm_dust %>% select(year, metric, value, dataset)
)

monthly_all <- bind_rows(
  monthly_all,
  monthly_pm_all  %>% select(month, metric, value, dataset),
  monthly_pm_dust %>% select(month, metric, value, dataset)
)

# -------------------- Plotting --------------------
metric_levels <- c("Wind speed (m/s)", "Visibility (km)",
                   "Relative humidity (%)", "PM10 (µg/m³)")
yearly_all  <- yearly_all  %>%
  mutate(metric  = fct_relevel(metric, metric_levels),
         dataset = factor(dataset, levels = c("Climatology","Dust events")))
monthly_all <- monthly_all %>%
  mutate(metric  = fct_relevel(metric, metric_levels),
         dataset = factor(dataset, levels = c("Climatology","Dust events")))

fill_cols   <- c(`Climatology` = "grey80", `Dust events` = "#2C7FB8")
stripe_cols <- c(`Climatology` = "grey80", `Dust events` = "#1B5E9A")

base_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  )

bar_plot <- function(d, x_var, ylab, xlab = NULL) {
  ggplot(
    d,
    aes(
      x = {{x_var}},
      y = value,
      fill = dataset,
      pattern = dataset,
      pattern_fill = dataset,
      pattern_angle = dataset   # map and scale below
    )
  ) +
    ggpattern::geom_col_pattern(
      position = position_dodge(width = 0.75),
      width = 0.65,
      colour = "grey20",
      pattern_colour = "grey20",
      pattern_density = 0.35,
      pattern_spacing = 0.02,
      pattern_size = 0.25
    ) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = fill_cols) +
    scale_pattern_manual(values = c(`Climatology` = "none", `Dust events` = "stripe")) +
    scale_pattern_fill_manual(values = stripe_cols) +
    scale_pattern_angle_manual(values = c(`Climatology` = 0, `Dust events` = 45)) +
    labs(y = ylab, x = xlab) +
    base_theme
}

# Yearly & Monthly panels
p_year <- bar_plot(yearly_all,  year,  ylab = "Yearly mean") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9))

p_mon  <- bar_plot(monthly_all, month, ylab = "Monthly mean", xlab = "Month")

final_plot <- p_year / p_mon
print(final_plot)

ggsave("wind_vis_rh_pm10_panels.png", final_plot, width = 18, height = 8, dpi = 300)
#---------------------------------------------------------------



# ---------------------------------------------------------------
# Packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(forcats)
  library(patchwork)
  library(ggpattern)   # install.packages("ggpattern")
  library(ggh4x)       # install.packages("ggh4x") for per-facet y limits
  library(readr)
})

# ---------------------------------------------------------------
# Station sets
stations_15 <- c("BFL","BAB","FAT","HJO","MAE","MYV","MCE","OVE",
                 "RBL","RDD","SAC","SMF","SCK","VIS","PTV")        # ASOS
pm_targets  <- c("HJO","BAB","BFL","FAT","MAE","MCE","MYV","SAC","SCK","VIS")  # PM10

# ---------------------------------------------------------------
# Helpers
prep_time <- function(df) {
  df %>%
    mutate(
      year  = year(local_time),
      month = factor(month.abb[month(local_time)],
                     levels = month.abb, ordered = TRUE),
      wind  = suppressWarnings(as.numeric(wind_speed_mps)),
      vis   = suppressWarnings(as.numeric(visibility_km)),
      rh    = suppressWarnings(as.numeric(RH_calc))
    )
}

yearly_summary <- function(df, label) {
  df %>%
    group_by(year) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-year, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

monthly_summary <- function(df, label) {
  df %>%
    group_by(month) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-month, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

mk_month_ord <- function(x) factor(as.character(x), levels = month.abb, ordered = TRUE)

# ---------------------------------------------------------------
# YOUR ASOS INPUTS MUST EXIST IN THE ENV:
#   dust_events_cv5 (dusty hours) and asos_data9 (all hours)
# ---------------------------------------------------------------
dust <- dust_events_cv5 %>%
  filter(station %in% stations_15) %>%
  prep_time()

yr_min <- min(dust$year, na.rm = TRUE)
yr_max <- max(dust$year, na.rm = TRUE)

climo <- asos_data9 %>%
  filter(station %in% stations_15,
         between(year(local_time), yr_min, yr_max)) %>%
  prep_time() %>%
  filter(is.na(wind) | wind > 0)   # remove calm hours only

# --- ASOS summaries (wind/vis/RH)
yearly_all  <- bind_rows(
  yearly_summary(climo, "Climatology"),
  yearly_summary(dust,  "Dust events")
)

monthly_all <- bind_rows(
  monthly_summary(climo, "Climatology"),
  monthly_summary(dust,  "Dust events")
)

# ---------------------------------------------------------------
# PM10 (add to the same panels)
# Requires pm_data.csv (Station, Date, PM10) and dusty-stations.csv (station, date [local])
# ---------------------------------------------------------------
pm_file   <- "pm_data.csv"
dust_file <- "dusty-stations.csv"

pm_all <- read_csv(pm_file, show_col_types = FALSE) %>%
  mutate(
    Station = toupper(trimws(Station)),
    dt      = parse_date_time(
      Date,
      orders = c("mdy HMS","mdy HM","mdy H","Ymd HMS","Y-m-d H:M:S","Y-m-d H:M","dmy HMS","dmy HM"),
      quiet  = TRUE
    ),
    date    = as_date(dt)
  ) %>%
  filter(Station %in% pm_targets, !is.na(PM10), !is.na(dt))

dust_keys <- dust_events_cv5 %>%
  transmute(
    Station = toupper(trimws(station)),
    date    = as_date(parse_date_time(date, orders = c("Y-m-d","Y/m/d","mdy"),
                                      exact = FALSE, quiet = TRUE))
  ) %>%
  filter(Station %in% pm_targets) %>%
  distinct()

pm_dust <- pm_all %>% inner_join(dust_keys, by = c("Station","date"))

# Yearly PM10
yearly_pm_all <- pm_all %>%
  mutate(year = year(dt)) %>%
  group_by(year) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Climatology")

yearly_pm_dust <- pm_dust %>%
  mutate(year = year(dt)) %>%
  group_by(year) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Dust events")

# Monthly PM10 (climatology across years)
monthly_pm_all <- pm_all %>%
  mutate(month = factor(month(dt), levels = 1:12, labels = month.abb, ordered = TRUE)) %>%
  group_by(month) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Climatology")

monthly_pm_dust <- pm_dust %>%
  mutate(month = factor(month(dt), levels = 1:12, labels = month.abb, ordered = TRUE)) %>%
  group_by(month) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Dust events")

# Harmonize month factor (fix bind_rows type mismatch)
monthly_pm_all  <- monthly_pm_all  %>% mutate(month = mk_month_ord(month))
monthly_pm_dust <- monthly_pm_dust %>% mutate(month = mk_month_ord(month))
monthly_all     <- monthly_all     %>% mutate(month = mk_month_ord(month))

# Append PM10 to the ASOS long tables
yearly_all  <- bind_rows(
  yearly_all,
  yearly_pm_all  %>% select(year, metric, value, dataset),
  yearly_pm_dust %>% select(year, metric, value, dataset)
)

monthly_all <- bind_rows(
  monthly_all,
  monthly_pm_all  %>% select(month, metric, value, dataset),
  monthly_pm_dust %>% select(month, metric, value, dataset)
)

# ---------------------------------------------------------------
# Plotting
metric_levels <- c("Wind speed (m/s)", "Visibility (km)",
                   "Relative humidity (%)", "PM10 (µg/m³)")

yearly_all  <- yearly_all  %>%
  mutate(metric  = fct_relevel(metric, metric_levels),
         dataset = factor(dataset, levels = c("Climatology","Dust events")))
monthly_all <- monthly_all %>%
  mutate(metric  = fct_relevel(metric, metric_levels),
         dataset = factor(dataset, levels = c("Climatology","Dust events")))

# Transparent red for Climatology; solid blue for Dust
fill_cols   <- c(`Climatology` = "#E53935", `Dust events` = "#2C7FB8")
stripe_cols <- c(`Climatology` = "#E53935", `Dust events` = "#1B5E9A")

base_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  )

bar_plot <- function(d, x_var, ylab, xlab = NULL, facet_y_scales = NULL) {
  p <- ggplot(
    d,
    aes(
      x = {{x_var}},
      y = value,
      fill = dataset,
      pattern = dataset,
      pattern_fill = dataset,
      alpha = dataset
    )
  ) +
    ggpattern::geom_col_pattern(
      position = position_dodge(width = 0.75),
      width = 0.65,
      colour = "grey20",
      pattern_colour = "grey20",
      pattern_density = 0.35,
      pattern_spacing = 0.02,
      pattern_size = 0.25
    ) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = fill_cols) +
    scale_alpha_manual(values = c(`Climatology` = 0.40, `Dust events` = 1.00),
                       guide = "none") +
    scale_pattern_manual(values = c(`Climatology` = "none", `Dust events` = "stripe")) +
    scale_pattern_fill_manual(values = stripe_cols) +
    labs(y = ylab, x = xlab) +
    base_theme
  
  if (!is.null(facet_y_scales)) {
    p <- p + facet_y_scales
  }
  p
}

# Yearly panel (free y)
p_year <- bar_plot(yearly_all, year, ylab = "Yearly mean") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9))

# Monthly panel with PM10 facet forced to 0–100
p_mon <- bar_plot(
  monthly_all, month, ylab = "Monthly mean", xlab = "Month",
  facet_y_scales = ggh4x::facetted_pos_scales(
    y = list(
      metric == "PM10 (µg/m³)" ~ scale_y_continuous(limits = c(0, 100))
    )
  )
)

final_plot <- p_year / p_mon
print(final_plot)

ggsave("wind_vis_rh_pm10_panels.png", final_plot, width = 18, height = 8, dpi = 300)
# ---------------------------------------------------------------




# ---------------------------------------------------------------
# Packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(forcats)
  library(patchwork)
  library(ggpattern)   # install.packages("ggpattern")
  library(ggh4x)       # install.packages("ggh4x")
  library(readr)
})

# ---------------------------------------------------------------
# Station sets
stations_15 <- c("BFL","BAB","FAT","HJO","MAE","MYV","MCE","OVE",
                 "RBL","RDD","SAC","SMF","SCK","VIS","PTV")          # ASOS
pm_targets  <- c("HJO","BAB","BFL","FAT","MAE","MCE","MYV","SAC","SCK","VIS")  # PM10

# ---------------------------------------------------------------
# Helpers
prep_time <- function(df) {
  df %>%
    mutate(
      year  = year(local_time),
      month = factor(month.abb[month(local_time)],
                     levels = month.abb, ordered = TRUE),
      wind  = suppressWarnings(as.numeric(wind_speed_mps)),
      vis   = suppressWarnings(as.numeric(visibility_km)),
      rh    = suppressWarnings(as.numeric(RH_calc))
    )
}

yearly_summary <- function(df, label) {
  df %>%
    group_by(year) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-year, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

monthly_summary <- function(df, label) {
  df %>%
    group_by(month) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-month, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

mk_month_ord <- function(x) factor(as.character(x), levels = month.abb, ordered = TRUE)

# ---------------------------------------------------------------
# ASOS input objects must exist: dust_events_cv5, asos_data9
# ---------------------------------------------------------------
dust <- dust_events_cv5 %>%
  filter(station %in% stations_15) %>%
  prep_time()

yr_min <- min(dust$year, na.rm = TRUE)
yr_max <- max(dust$year, na.rm = TRUE)

climo <- asos_data9 %>%
  filter(station %in% stations_15,
         between(year(local_time), yr_min, yr_max)) %>%
  prep_time() %>%
  filter(is.na(wind) | wind > 0)   # remove calm hours only

# ASOS summaries
yearly_all  <- bind_rows(
  yearly_summary(climo, "Climatology"),
  yearly_summary(dust,  "Dust events")
)

monthly_all <- bind_rows(
  monthly_summary(climo, "Climatology"),
  monthly_summary(dust,  "Dust events")
)

# ---------------------------------------------------------------
# PM10 (join with dusty days and add to panels)
# Requires: pm_data.csv (Station, Date, PM10) and dusty-stations.csv (station, date [local])
# ---------------------------------------------------------------
pm_file   <- "pm_data.csv"
dust_file <- "dusty-stations.csv"

pm_all <- read_csv(pm_file, show_col_types = FALSE) %>%
  mutate(
    Station = toupper(trimws(Station)),
    dt      = parse_date_time(
      Date,
      orders = c("mdy HMS","mdy HM","mdy H","Ymd HMS","Y-m-d H:M:S","Y-m-d H:M","dmy HMS","dmy HM"),
      quiet  = TRUE
    ),
    date    = as_date(dt)
  ) %>%
  filter(Station %in% pm_targets, !is.na(PM10), !is.na(dt))

dust_keys <- dust_events_cv5 %>%
  transmute(
    Station = toupper(trimws(station)),
    date    = as_date(parse_date_time(date, orders = c("Y-m-d","Y/m/d","mdy"),
                                      exact = FALSE, quiet = TRUE))
  ) %>%
  filter(Station %in% pm_targets) %>%
  distinct()

pm_dust <- pm_all %>% inner_join(dust_keys, by = c("Station","date"))

# Yearly PM10
yearly_pm_all <- pm_all %>%
  mutate(year = year(dt)) %>%
  group_by(year) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Climatology")

yearly_pm_dust <- pm_dust %>%
  mutate(year = year(dt)) %>%
  group_by(year) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Dust events")

# Monthly PM10
monthly_pm_all <- pm_all %>%
  mutate(month = factor(month(dt), levels = 1:12, labels = month.abb, ordered = TRUE)) %>%
  group_by(month) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Climatology")

monthly_pm_dust <- pm_dust %>%
  mutate(month = factor(month(dt), levels = 1:12, labels = month.abb, ordered = TRUE)) %>%
  group_by(month) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Dust events")

# Harmonize month factor (fix bind_rows mismatch)
monthly_pm_all  <- monthly_pm_all  %>% mutate(month = mk_month_ord(month))
monthly_pm_dust <- monthly_pm_dust %>% mutate(month = mk_month_ord(month))
monthly_all     <- monthly_all     %>% mutate(month = mk_month_ord(month))

# Append PM10 into long tables
yearly_all  <- bind_rows(
  yearly_all,
  yearly_pm_all  %>% select(year, metric, value, dataset),
  yearly_pm_dust %>% select(year, metric, value, dataset)
)

monthly_all <- bind_rows(
  monthly_all,
  monthly_pm_all  %>% select(month, metric, value, dataset),
  monthly_pm_dust %>% select(month, metric, value, dataset)
)

# ---------------------------------------------------------------
# Plotting
metric_levels <- c("Wind speed (m/s)", "Visibility (km)",
                   "Relative humidity (%)", "PM10 (µg/m³)")

yearly_all  <- yearly_all  %>%
  mutate(metric  = fct_relevel(metric, metric_levels),
         dataset = factor(dataset, levels = c("Climatology","Dust events")))
monthly_all <- monthly_all %>%
  mutate(metric  = fct_relevel(metric, metric_levels),
         dataset = factor(dataset, levels = c("Climatology","Dust events")))

# Outline-only red for Climatology; solid blue (striped) for Dust
fill_vals   <- c(`Climatology` = NA,          `Dust events` = "#2C7FB8")
border_vals <- c(`Climatology` = "#E53935",   `Dust events` = "grey20")
pattern_vals      <- c(`Climatology` = "none",  `Dust events` = "stripe")
pattern_fill_vals <- c(`Climatology` = NA,      `Dust events` = "#1B5E9A")

base_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  )

bar_plot <- function(d, x_var, ylab, xlab = NULL, facet_y_scales = NULL) {
  p <- ggplot(
    d,
    aes(
      x = {{x_var}},
      y = value,
      group = dataset,
      fill = dataset,
      colour = dataset,        # outline mapping
      pattern = dataset,
      pattern_fill = dataset
    )
  ) +
    ggpattern::geom_col_pattern(
      position = position_dodge(width = 0.75),
      width = 0.65,
      # outline aesthetics come from mapped colour
      pattern_colour = "grey20",
      pattern_density = 0.35,
      pattern_spacing = 0.02,
      pattern_size = 0.25
    ) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    # NOTE: NA fill for Climatology gives outline-only bars
    scale_fill_manual(values = fill_vals, na.translate = FALSE) +
    scale_colour_manual(values = border_vals, guide = "none") +
    scale_pattern_manual(values = pattern_vals) +
    scale_pattern_fill_manual(values = pattern_fill_vals) +
    labs(y = ylab, x = xlab) +
    base_theme
  
  # tweak legend keys so Climatology shows a red outline
  p <- p + guides(
    fill = guide_legend(override.aes = list(
      fill = c(NA, "#2C7FB8"),
      colour = c("#E53935", "grey20"),
      pattern = c("none","stripe")
    ))
  )
  
  if (!is.null(facet_y_scales)) {
    p <- p + facet_y_scales
  }
  p
}

# Yearly panel
p_year <- bar_plot(yearly_all, year, ylab = "Yearly mean") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9))

# Monthly panel with PM10 facet forced to 0–100
p_mon <- bar_plot(
  monthly_all, month, ylab = "Monthly mean", xlab = "Month",
  facet_y_scales = ggh4x::facetted_pos_scales(
    y = list(metric == "PM10 (µg/m³)" ~ scale_y_continuous(limits = c(0, 100)))
  )
)

final_plot <- p_year / p_mon
print(final_plot)

ggsave("wind_vis_rh_pm10_panels.png", final_plot, width = 18, height = 8, dpi = 500)
# ---------------------------------------------------------------



# ---------------------------------------------------------------
# Full script: Wind / Visibility / RH + PM10
# - Climatology vs Dust events
# - BOTH series transparent (outline-only); Dust has blue stripes
# - PM10 monthly facet fixed to 0–100
# - Fixes month factor mismatch for bind_rows()
#
# Requirements in working dir:
#   pm_data.csv (columns: Station, Date, PM10)
#   dusty-stations.csv (columns: station, date [LOCAL time])
#
# Assumes in memory:
#   dust_events_cv5 (ASOS rows for dusty hours)
#   asos_data9      (ASOS rows for all hours)
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(forcats)
  library(patchwork)
  library(ggpattern)   # install.packages("ggpattern")
  library(ggh4x)       # install.packages("ggh4x")
  library(readr)
})

# -------------------- Station sets --------------------
stations_15 <- c("BFL","BAB","FAT","HJO","MAE","MYV","MCE","OVE",
                 "RBL","RDD","SAC","SMF","SCK","VIS","PTV")     # ASOS vars
pm_targets  <- c("HJO","BAB","BFL","FAT","MAE","MCE","MYV","SAC","SCK","VIS")  # PM10

# -------------------- Helpers --------------------
prep_time <- function(df) {
  df %>%
    mutate(
      year  = year(local_time),
      month = factor(month.abb[month(local_time)], levels = month.abb, ordered = TRUE),
      wind  = suppressWarnings(as.numeric(wind_speed_mps)),
      vis   = suppressWarnings(as.numeric(visibility_km)),
      rh    = suppressWarnings(as.numeric(RH_calc))
    )
}

yearly_summary <- function(df, label) {
  df %>%
    group_by(year) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-year, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

monthly_summary <- function(df, label) {
  df %>%
    group_by(month) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-month, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

mk_month_ord <- function(x) factor(as.character(x), levels = month.abb, ordered = TRUE)

# -------------------- Build ASOS datasets --------------------
dust <- dust_events_cv5 %>%
  filter(station %in% stations_15) %>%
  prep_time()

yr_min <- min(dust$year, na.rm = TRUE)
yr_max <- max(dust$year, na.rm = TRUE)

climo <- asos_data9 %>%
  filter(station %in% stations_15,
         between(year(local_time), yr_min, yr_max)) %>%
  prep_time() %>%
  filter(is.na(wind) | wind > 0)   # remove calm hours only

yearly_all  <- bind_rows(
  yearly_summary(climo, "Climatology"),
  yearly_summary(dust,  "Dust events")
)

monthly_all <- bind_rows(
  monthly_summary(climo, "Climatology"),
  monthly_summary(dust,  "Dust events")
)

# -------------------- PM10: read, join to dusty days, summarise --------------------
pm_file   <- "pm_data.csv"
dust_file <- "dusty-stations.csv"

pm_all <- read_csv(pm_file, show_col_types = FALSE) %>%
  mutate(
    Station = toupper(trimws(Station)),
    dt      = parse_date_time(
      Date,
      orders = c("mdy HMS","mdy HM","mdy H","Ymd HMS","Y-m-d H:M:S","Y-m-d H:M","dmy HMS","dmy HM"),
      quiet  = TRUE
    ),
    date    = as_date(dt)
  ) %>%
  filter(Station %in% pm_targets, !is.na(PM10), !is.na(dt))

dust_keys <- dust_events_cv5%>%
  transmute(
    Station = toupper(trimws(station)),
    date    = as_date(parse_date_time(date, orders = c("Y-m-d","Y/m/d","mdy"),
                                      exact = FALSE, quiet = TRUE))
  ) %>%
  filter(Station %in% pm_targets) %>%
  distinct()

pm_dust <- pm_all %>% inner_join(dust_keys, by = c("Station","date"))

# Yearly PM10
yearly_pm_all <- pm_all %>%
  mutate(year = year(dt)) %>%
  group_by(year) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Climatology")

yearly_pm_dust <- pm_dust %>%
  mutate(year = year(dt)) %>%
  group_by(year) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Dust events")

# Monthly PM10
monthly_pm_all <- pm_all %>%
  mutate(month = factor(month(dt), levels = 1:12, labels = month.abb, ordered = TRUE)) %>%
  group_by(month) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Climatology")

monthly_pm_dust <- pm_dust %>%
  mutate(month = factor(month(dt), levels = 1:12, labels = month.abb, ordered = TRUE)) %>%
  group_by(month) %>%
  summarise(value = mean(PM10, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = "PM10 (µg/m³)", dataset = "Dust events")

# Harmonize month factor for bind_rows
monthly_pm_all  <- monthly_pm_all  %>% mutate(month = mk_month_ord(month))
monthly_pm_dust <- monthly_pm_dust %>% mutate(month = mk_month_ord(month))
monthly_all     <- monthly_all     %>% mutate(month = mk_month_ord(month))

# Append PM10 into long tables
yearly_all  <- bind_rows(
  yearly_all,
  yearly_pm_all  %>% select(year, metric, value, dataset),
  yearly_pm_dust %>% select(year, metric, value, dataset)
)

monthly_all <- bind_rows(
  monthly_all,
  monthly_pm_all  %>% select(month, metric, value, dataset),
  monthly_pm_dust %>% select(month, metric, value, dataset)
)

# -------------------- Plotting (transparent outlines + dust stripes) --------------------
metric_levels <- c("Wind speed (m/s)", "Visibility (km)",
                   "Relative humidity (%)", "PM10 (µg/m³)")

yearly_all  <- yearly_all  %>%
  mutate(metric  = fct_relevel(metric, metric_levels),
         dataset = factor(dataset, levels = c("Climatology","Dust events")))
monthly_all <- monthly_all %>%
  mutate(metric  = fct_relevel(metric, metric_levels),
         dataset = factor(dataset, levels = c("Climatology","Dust events")))

# Transparent (no fill) for both; outlines differ; dust has stripes
fill_vals           <- c(`Climatology` = NA,         `Dust events` = NA)          # transparent
border_vals         <- c(`Climatology` = "#E53935",  `Dust events` = "#1B5E9A")   # red vs blue outlines
pattern_vals        <- c(`Climatology` = "none",     `Dust events` = "stripe")
pattern_fill_vals   <- c(`Climatology` = NA,         `Dust events` = NA)          # transparent pattern bg
pattern_colour_vals <- c(`Climatology` = "#E53935",  `Dust events` = "#1B5E9A")   # stripe color

base_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  )

bar_plot <- function(d, x_var, ylab, xlab = NULL, facet_y_scales = NULL) {
  p <- ggplot(
    d,
    aes(
      x = {{x_var}},
      y = value,
      group = dataset,
      fill = dataset,             # both NA (transparent)
      colour = dataset,           # outlines
      pattern = dataset,          # stripes only for Dust events
      pattern_fill = dataset,     # transparent
      pattern_colour = dataset    # stripe color
    )
  ) +
    ggpattern::geom_col_pattern(
      position = position_dodge(width = 0.75),
      width = 0.65,
      pattern_density = 0.35,
      pattern_spacing = 0.02,
      pattern_size    = 0.25
    ) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = fill_vals, na.translate = FALSE) +
    scale_colour_manual(values = border_vals, guide = "none") +
    scale_pattern_manual(values = pattern_vals) +
    scale_pattern_fill_manual(values = pattern_fill_vals) +
    scale_pattern_colour_manual(values = pattern_colour_vals, guide = "none") +
    labs(y = ylab, x = xlab) +
    base_theme +
    guides(
      fill = guide_legend(override.aes = list(
        fill          = c(NA, NA),
        colour        = c("#E53935", "#1B5E9A"),
        pattern       = c("none","stripe"),
        pattern_fill  = c(NA, NA),
        pattern_colour= c("#E53935", "#1B5E9A")
      ))
    )
  
  if (!is.null(facet_y_scales)) {
    p <- p + facet_y_scales
  }
  p
}

# Yearly panel
p_year <- bar_plot(yearly_all, year, ylab = "Yearly mean") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9))

# Monthly panel with PM10 facet forced to 0–100
p_mon <- bar_plot(
  monthly_all, month, ylab = "Monthly mean", xlab = "Month",
  facet_y_scales = ggh4x::facetted_pos_scales(
    y = list(metric == "PM10 (µg/m³)" ~ scale_y_continuous(limits = c(0, 100)))
  )
)

final_plot <- p_year / p_mon
print(final_plot)

ggsave("wind_vis_rh_pm10_panels.png", final_plot, width = 18, height = 8, dpi = 500)
# ---------------------------------------------------------------










#--------------------------------------------------------------------
# ===========================
#  PDFS FOR DUST vs CLIMO
# ===========================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(forcats)
  library(patchwork)   # to arrange panels
})

# ----- Inputs expected in memory -----
# dust_events_cv5: station, local_time, wind_speed_mps, visibility_km, RH_calc
# asos_datat9 (or asos_data9): same fields for climatology

# Pick whichever climatology object you actually have
climo_raw <- if (exists("asos_datat9")) {
  asos_datat9
} else if (exists("asos_data9")) {
  asos_data9
} else {
  stop("Please provide `asos_datat9` or `asos_data9` in the workspace.")
}

# ------------------ Helpers ------------------
prep_time <- function(df) {
  df %>%
    transmute(
      station,
      local_time = as.POSIXct(local_time),
      date  = as.Date(local_time),
      month = factor(month.abb[month(local_time)], levels = month.abb, ordered = TRUE),
      wind  = suppressWarnings(as.numeric(wind_speed_mps)),
      vis   = suppressWarnings(as.numeric(visibility_km)),
      rh    = suppressWarnings(as.numeric(RH_calc))
    )
}

# Harmonize visibility: cap at common upper bound so ">=10 km" / ">=10 sm" are comparable
harmonize_vis <- function(x, cap_km = 16) {
  x <- as.numeric(x)
  x <- ifelse(is.na(x), NA_real_, ifelse(x > cap_km, cap_km, x))
  x
}

# Collapse to one value per station-day so each day contributes once
collapse_station_day <- function(df) {
  df %>%
    group_by(station, date, month) %>%
    summarise(
      wind = median(wind, na.rm = TRUE),
      vis  = median(vis,  na.rm = TRUE),
      rh   = median(rh,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(c(wind, vis, rh), names_to = "metric_raw", values_to = "value") %>%
    mutate(metric = recode(metric_raw,
                           wind = "Wind speed (m/s)",
                           vis  = "Visibility (km)",
                           rh   = "Relative humidity (%)"))
}

# ------------------ Prep data ------------------
dust  <- dust_events_cv5 %>% prep_time()

# Restrict climatology to the dust time span and overlapping stations
yrs  <- range(year(dust$local_time), na.rm = TRUE)
stns <- intersect(unique(dust$station), unique(climo_raw$station))

climo <- climo_raw %>%
  prep_time() %>%
  filter(station %in% stns,
         between(year(local_time), yrs[1], yrs[2]),
         is.na(wind) | wind > 0)              # remove calm hours (keep NA)

# Harmonize visibility & basic sanity limits
dust$vis  <- harmonize_vis(dust$vis)
climo$vis <- harmonize_vis(climo$vis)

# Keep plausible ranges
clip_ok <- function(d) {
  d %>%
    mutate(
      wind = ifelse(wind >= 0 & is.finite(wind), wind, NA_real_),
      vis  = ifelse(vis  >= 0 & is.finite(vis),  vis,  NA_real_),
      rh   = ifelse(rh   >= 0 & rh <= 100,       rh,   NA_real_)
    )
}
dust  <- clip_ok(dust)
climo <- clip_ok(climo)

# Collapse to station-day
dust_d  <- collapse_station_day(dust)  %>% mutate(dataset = "Dust events")
climo_d <- collapse_station_day(climo) %>% mutate(dataset = "Climatology")

pdf_daily <- bind_rows(dust_d, climo_d) %>%
  filter(is.finite(value)) %>%
  mutate(metric = fct_relevel(metric,
                              c("Wind speed (m/s)", "Visibility (km)", "Relative humidity (%)")))

# ------------------ Plotting helpers ------------------
fill_cols <- c(`Climatology` = "grey80", `Dust events` = "#2C7FB8")
line_cols <- c(`Climatology` = "grey35", `Dust events` = "#1B5E9A")

plot_overall <- function(m) {
  d <- pdf_daily %>% filter(metric == m)
  xr <- quantile(d$value, c(.01, .99), na.rm = TRUE)
  meds <- d %>% group_by(dataset) %>% summarise(med = median(value, na.rm = TRUE), .groups = "drop")
  ggplot(d, aes(value, fill = dataset, colour = dataset)) +
    geom_density(alpha = .28, adjust = 2.5, linewidth = .8) +   # wider bandwidth for smoothness
    geom_vline(data = meds, aes(xintercept = med, colour = dataset),
               linetype = 2, linewidth = .8, show.legend = FALSE) +
    coord_cartesian(xlim = xr) +
    scale_fill_manual(values = fill_cols) +
    scale_colour_manual(values = line_cols) +
    labs(x = m, y = "Density", title = m) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "top",
          legend.title    = element_blank())
}

plot_overall <- function(m) {
  d <- pdf_daily %>% dplyr::filter(metric == m)
  
  # fixed axis ranges + nice ticks per metric
  ax <- switch(
    m,
    "Wind speed (m/s)"       = list(lim = c(0, 12),  brk = seq(0, 12, 2)),
    "Visibility (km)"        = list(lim = c(0, 20),  brk = seq(0, 20, 5)),
    "Relative humidity (%)"  = list(lim = c(0, 100), brk = seq(0, 100, 25)),
    # fallback if a different metric name sneaks in
    {rng <- quantile(d$value, c(.01, .99), na.rm = TRUE)
    list(lim = rng, brk = waiver())}
  )
  
  meds <- d %>%
    dplyr::group_by(dataset) %>%
    dplyr::summarise(med = median(value, na.rm = TRUE), .groups = "drop")
  
  ggplot(d, aes(x = value, fill = dataset, colour = dataset)) +
    geom_density(alpha = .28, adjust = 2.5, linewidth = .8) +
    geom_vline(data = meds, aes(xintercept = med, colour = dataset),
               linetype = 2, linewidth = .8, show.legend = FALSE) +
    coord_cartesian(xlim = ax$lim) +
    scale_x_continuous(breaks = ax$brk, expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = fill_cols, labels = c("Climatology","Dust events")) +
    scale_colour_manual(values = line_cols, guide = "none") +
    labs(x = m, y = "Density") +             # no title
    theme_minimal(base_size = 19) +
    theme(
      legend.position = "top",
      legend.title    = element_blank(),
      plot.title      = element_blank()
    )
}


plot_monthly <- function(m, ncol = 4) {
  d <- pdf_daily %>% filter(metric == m)
  xr <- quantile(d$value, c(.01, .99), na.rm = TRUE)
  
  ggplot(d, aes(value, fill = dataset, colour = dataset)) +
    geom_density(alpha = .28, adjust = 2.2, linewidth = .6) +
    coord_cartesian(xlim = xr) +
    facet_wrap(~ month, ncol = ncol, scales = "free_y") +   # << key change
    scale_fill_manual(values = fill_cols) +
    scale_colour_manual(values = line_cols) +
    labs(x = m, y = "Density") +
    theme_minimal(base_size = 18) +
    theme(
      legend.position = "top",
      legend.title    = element_blank(),
      panel.grid.minor = element_blank()
    )
}



# ------------------ Draw & (optionally) save ------------------
p_overall <- plot_overall("Wind speed (m/s)") |
  plot_overall("Visibility (km)")  |
  plot_overall("Relative humidity (%)")
p_overall

# ------------------ Draw & (optionally) save ------------------
p_overall <- plot_overall("Wind speed (m/s)") |
  plot_overall("Visibility (km)")  |
  plot_overall("Relative humidity (%)")
p_overall

ggsave("pdf_overall_dust_vs_climo.png", p_overall, width = 24, height = 10, dpi = 300)

p_wind_m  <- plot_monthly("Wind speed (m/s)")
p_vis_m   <- plot_monthly("Visibility (km)")
p_rh_m    <- plot_monthly("Relative humidity (%)")

p_wind_m  # not cramped: 4 columns
p_vis_m
p_rh_m

plot_overall2 <- p_wind_m / p_vis_m / p_rh_m
plot_overall2 <- plot_overall2 & theme(plot.title = element_blank())


 
ggsave("pdf_monthly_wind.png", p_wind_m, width = 14, height = 8, dpi = 300)
ggsave("pdf_monthly_vis.png",  p_vis_m,  width = 14, height = 8, dpi = 300)
ggsave("pdf_monthly_rh.png",   p_rh_m,   width = 14, height = 8, dpi = 300)


final_plot2 <- p_mon/p_overall 

final_plot2
ggsave("combined_pattern.png", final_plot2, width = 18, height = 8, dpi = 300)

ggsave("mon1_pattern.png", p_mon, width = 18, height = 8, dpi = 300)
ggsave("mon2_pattern.png", plot_overall2, width = 18, height = 8, dpi = 300)

ggsave("mon2_pattern.png", plot_overall2, width = 18, height = 21, dpi = 300)









#--------------precipitation---------------------
# --- Helpers ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(forcats)
library(ggplot2)

# Colors to match your other plots
fill_cols <- c(`Climatology` = "grey80", `Dust events` = "#2C7FB8")
line_cols <- c(`Climatology` = "grey30", `Dust events` = "#1B5E9A")

# Clean precip column like p01m (strings "null", blanks, "T" -> 0 or NA)
clean_precip <- function(x) {
  x <- as.character(x)
  x[x %in% c("null","NULL","","NaN","NA")] <- NA
  # Treat trace as 0 (change to 0.005 if you prefer)
  x[x %in% c("T","t","TRACE","Trace","trace")] <- "0"
  suppressWarnings(as.numeric(x))
}

# Pick a precipitation column if names differ across files
pick_precip <- function(df) {
  if ("p01m" %in% names(df)) clean_precip(df$p01m)
  else if ("precip_mm" %in% names(df)) clean_precip(df$precip_mm)
  else if ("precip" %in% names(df)) clean_precip(df$precip)
  else rep(NA_real_, nrow(df))
}

# Add time-derived fields + variables (extends what you already do)
prep_time_with_prcp <- function(df) {
  df %>%
    mutate(
      date  = as.Date(local_time),
      year  = year(local_time),
      month = factor(month.abb[month(local_time)], levels = month.abb, ordered = TRUE),
      wind  = suppressWarnings(as.numeric(wind_speed_mps)),
      vis   = suppressWarnings(as.numeric(visibility_km)),
      rh    = suppressWarnings(as.numeric(RH_calc)),
      prcp  = pick_precip(cur_data())
      # If your precip is in inches, convert to mm:
      # prcp  = pick_precip(cur_data()) * 25.4
    )
}

# --- Data prep (reuses your existing objects) ------------------------------
# 15 stations subset (adjust to your set)
stations_15 <- c("BFL","BAB","FAT","HJO","MAE","MYV","MCE",
                 "OVE","RBL","RDD","SAC","SMF","SCK","VIS","PTV")

# Dust events
dust <- dust_events_cv7 %>%
  filter(station %in% stations_15) %>%
  prep_time_with_prcp()

# Match climatology period to dust years, filter calm hours as before
yr_min <- min(dust$year, na.rm = TRUE)
yr_max <- max(dust$year, na.rm = TRUE)

climo <- asos_data19 %>%
  filter(station %in% stations_15,
         between(year(local_time), yr_min, yr_max)) %>%
  prep_time_with_prcp() %>%
  # keep NA, remove only *calm* hours for wind as you requested earlier
  filter(is.na(wind) | wind > 0)

# Build a tidy frame for precipitation PDFs
pdf_prcp <- bind_rows(
  dust  %>% mutate(dataset = "Dust events"),
  climo %>% mutate(dataset = "Climatology")
) %>%
  select(date, month, dataset, prcp) %>%
  rename(value = prcp) %>%
  mutate(metric = "Precipitation (mm)")              # change label/unit if needed

# Handy check: proportion of dry hours (zero precip) by dataset & month
dry_share_overall <- pdf_prcp %>%
  summarise(dry_share = mean(value <= 0 | is.na(value), na.rm = TRUE),
            .by = dataset)

dry_share_monthly <- pdf_prcp %>%
  summarise(dry_share = mean(value <= 0 | is.na(value), na.rm = TRUE),
            .by = c(month, dataset))

# --- Plot functions (overall & monthly) ------------------------------------
# Overall precipitation PDF (positive precipitation only)
plot_overall_precip <- function() {
  d_all <- pdf_prcp %>% filter(!is.na(value) & value > 0)
  if (nrow(d_all) == 0) {
    stop("No positive precipitation found after cleaning.")
  }
  xr <- c(0, quantile(d_all$value, 0.995, na.rm = TRUE))
  meds <- d_all %>%
    summarise(med = median(value, na.rm = TRUE), .by = dataset)
  
  ggplot(d_all, aes(value, fill = dataset, colour = dataset)) +
    geom_density(alpha = 0.28, adjust = 2.0, linewidth = 0.8) +
    geom_vline(data = meds, aes(xintercept = med, colour = dataset),
               linetype = 2, linewidth = 0.9, show.legend = FALSE) +
    coord_cartesian(xlim = xr) +
    scale_fill_manual(values = fill_cols) +
    scale_colour_manual(values = line_cols) +
    labs(x = "Precipitation (mm)", y = "Density") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    )
}

# Monthly precipitation PDFs (positive precipitation only, free y)
plot_monthly_precip <- function(ncol = 4) {
  d_all <- pdf_prcp %>% filter(!is.na(value) & value > 0)
  if (nrow(d_all) == 0) {
    stop("No positive precipitation found after cleaning.")
  }
  ggplot(d_all, aes(value, fill = dataset, colour = dataset)) +
    geom_density(alpha = 0.25, adjust = 1.8, linewidth = 0.6) +
    scale_fill_manual(values = fill_cols) +
    scale_colour_manual(values = line_cols) +
    facet_wrap(~ month, ncol = ncol, scales = "free_y") +
    labs(x = "Precipitation (mm)", y = "Density") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# --- Make the figures ------------------------------------------------------
p_precip_overall <- plot_overall_precip()
p_precip_monthly <- plot_monthly_precip()

p_precip_overall
p_precip_monthly

# If you want, annotate dry-hour share next to/under the plots:
dry_share_overall %>%
  mutate(`Dry or zero precip (%)` = round(100 * dry_share, 1)) %>%
  select(dataset, `Dry or zero precip (%)`) %>%
  print()







suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(fitdistrplus); library(scales)
})

# expects: pdf_prcp with columns: dataset (factor: "Climatology","Dust events"), value (mm, >=0)
# and palettes you already used:
fill_cols <- c(`Climatology` = "grey80", `Dust events` = "#9ecae1")
line_cols <- c(`Climatology` = "grey30", `Dust events` = "#2c7fb8")

fit_zig <- function(x) {
  x <- x[!is.na(x) & x >= 0]
  p0 <- mean(x == 0)
  xp <- x[x > 0]
  if (length(xp) >= 2 && sum(xp) > 0) {
    fit <- suppressWarnings(fitdistrplus::fitdist(xp, "gamma"))
    shape <- unname(fit$estimate["shape"])
    rate  <- unname(fit$estimate["rate"])
  } else {
    # fallback if almost all zeros: a tiny exponential-ish tail
    shape <- 1
    rate  <- 1 / pmax(mean(xp), 0.1)
  }
  list(p0 = p0, shape = shape, rate = rate)
}

zig_median <- function(p0, shape, rate) {
  if (p0 >= 0.5) 0 else qgamma((0.5 - p0) / (1 - p0), shape = shape, rate = rate)
}

plot_precip_pdf_parametric <- function(pdf_prcp) {
  d <- pdf_prcp %>% filter(!is.na(value), value >= 0)
  
  # fit separate ZIG for each dataset
  fits <- d %>%
    group_by(dataset) %>%
    summarise(mod = list(fit_zig(value)), .groups = "drop") %>%
    mutate(
      p0    = purrr::map_dbl(mod, "p0"),
      shape = purrr::map_dbl(mod, "shape"),
      rate  = purrr::map_dbl(mod, "rate"),
      med   = mapply(zig_median, p0, shape, rate),
      .keep = "unused"
    )
  
  # x-range: 0 to 99.5th percentile of positive amounts in the combined data
  pos_all <- d$value[d$value > 0]
  x_hi <- if (length(pos_all)) quantile(pos_all, 0.995, na.rm = TRUE) else 1
  x_grid <- seq(0, x_hi, length.out = 1000)
  
  # parametric continuous PDF (Gamma part scaled by 1-p0)
  curves <- fits %>%
    tidyr::expand_grid(x = x_grid) %>%
    mutate(pdf = (1 - p0) * dgamma(x, shape = shape, rate = rate))
  
  # choose a spike height (purely visual—Dirac mass has infinite height)
  spike_heights <- curves %>%
    group_by(dataset) %>%
    summarise(h = max(pdf) * 1.05, .groups = "drop") %>%
    left_join(fits %>% select(dataset, p0), by = "dataset")
  
  ggplot() +
    # light histogram to give scale (no kernel density)
    geom_histogram(
      data = d,
      aes(x = value, y = after_stat(density), fill = dataset),
      bins = 50, position = "identity", alpha = 0.18, colour = NA
    ) +
    # parametric PDF (Gamma part)
    geom_area(
      data = curves,
      aes(x = x, y = pdf, fill = dataset),
      alpha = 0.10, colour = NA
    ) +
    geom_line(
      data = curves,
      aes(x = x, y = pdf, colour = dataset),
      linewidth = 1.1
    ) +
    # visual spike at zero + label P(0)
    geom_segment(
      data = spike_heights,
      aes(x = 0, xend = 0, y = 0, yend = h, colour = dataset),
      linewidth = 1.2
    ) +
    geom_text(
      data = spike_heights,
      aes(x = 0, y = h, label = paste0("P(0)=", scales::percent(p0, accuracy = 0.1)),
          colour = dataset),
      hjust = -0.05, vjust = 1.1, size = 3.8, show.legend = FALSE
    ) +
    # mixture median (vertical, overall median of ZIG)
    geom_vline(
      data = fits, aes(xintercept = med, colour = dataset),
      linetype = 2, linewidth = 0.9, show.legend = FALSE
    ) +
    coord_cartesian(xlim = c(0, x_hi)) +
    scale_fill_manual(values = fill_cols) +
    scale_colour_manual(values = line_cols) +
    labs(x = "Precipitation (mm)", y = "Density") +
    theme_minimal(base_size = 24) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}


p_prcp <- plot_precip_pdf_parametric(pdf_prcp)
p_prcp

ggsave("precipitation5.png", p_prcp, width = 18, height = 9, dpi = 300)









#----------------------------------------------

# -------------------- Packages --------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(forcats)
library(patchwork)

# -------------------- Config ----------------------
stations_15 <- c("BFL","BAB","FAT","HJO","MAE","MYV","MCE","OVE",
                 "RBL","RDD","SAC","SMF","SCK","VIS","PTV")

# Helper: add time fields and coerce numerics
prep_time <- function(df) {
  df %>%
    mutate(
      year  = year(local_time),
      month = factor(month.abb[month(local_time)],
                     levels = month.abb, ordered = TRUE),
      wind  = suppressWarnings(as.numeric(wind_speed_mps)),
      vis   = suppressWarnings(as.numeric(visibility_km)),
      rh    = suppressWarnings(as.numeric(RH_calc))
    )
}

# --- Precipitation helpers (keep zeros) -------------------------------
.pick_first <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  rep(NA, nrow(df))
}
.parse_prcp <- function(x) {
  v <- as.character(x)
  v[is.na(v) | v == "" | grepl("^null$", v, ignore.case = TRUE)] <- "0"
  suppressWarnings(as.numeric(gsub("[^0-9.]+", "", v)))
}
add_precip <- function(df) {
  cand <- c("p01m","p01h","precip_mm","precip","prcp_mm","prcp","P01M","P01I")
  df$prcp <- .parse_prcp(.pick_first(df, cand))
  df
}

# -------------------- Data prep --------------------
# Dust dataset (already dust-only rows)
dust <- dust_events_cv7 %>%
  filter(station %in% stations_15) %>%
  prep_time() %>%
  add_precip()

# Year window from dust data
yr_min <- min(dust$year, na.rm = TRUE)
yr_max <- max(dust$year, na.rm = TRUE)

# Climatology: same stations/years, remove calm (0 m/s) but keep NA wind
climo <- asos_data19 %>%
  filter(station %in% stations_15,
         between(year(local_time), yr_min, yr_max)) %>%
  prep_time() %>%
  filter(is.na(wind) | wind > 0) %>%
  add_precip()

# -------------------- Summaries (now including precip) ----------------
yearly_summary <- function(df, label) {
  df %>%
    group_by(year) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      `Precipitation (mm)`    = mean(prcp, na.rm = TRUE),   # zeros kept
      .groups = "drop"
    ) %>%
    pivot_longer(-year, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

monthly_summary <- function(df, label) {
  df %>%
    group_by(month) %>%
    summarise(
      `Wind speed (m/s)`      = mean(wind, na.rm = TRUE),
      `Visibility (km)`       = mean(vis,  na.rm = TRUE),
      `Relative humidity (%)` = mean(rh,   na.rm = TRUE),
      `Precipitation (mm)`    = mean(prcp, na.rm = TRUE),   # zeros kept
      .groups = "drop"
    ) %>%
    pivot_longer(-month, names_to = "metric", values_to = "value") %>%
    mutate(dataset = label)
}

yearly_all  <- bind_rows(
  yearly_summary(climo, "Climatology"),
  yearly_summary(dust,  "Dust events")
)

monthly_all <- bind_rows(
  monthly_summary(climo, "Climatology"),
  monthly_summary(dust,  "Dust events")
)

# Optional: control facet order
metric_levels <- c("Wind speed (m/s)", "Visibility (km)",
                   "Relative humidity (%)", "Precipitation (mm)")
yearly_all  <- yearly_all  %>% mutate(metric = fct_relevel(metric, metric_levels))
monthly_all <- monthly_all %>% mutate(metric = fct_relevel(metric, metric_levels))

# -------------------- Plotting ---------------------
pal <- c(`Climatology` = "grey80", `Dust events` = "#2C7FB8")

base_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  )

# Yearly (top row)
p_year <- ggplot(yearly_all, aes(x = year, y = value, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pal) +
  labs(y = "Yearly mean") +
  base_theme +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))

# Monthly (bottom row)
p_mon <- ggplot(monthly_all, aes(x = month, y = value, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pal) +
  labs(y = "Monthly mean", x = "Month") +
  base_theme

final_plot <- p_year / p_mon
final_plot

# Save if you like
ggsave("climo_vs_dust_year_month_with_precip.png", final_plot, width = 14, height = 8, dpi = 300)






#---------------------------------------------





# ---- Packages ----
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lubridate)
  library(ggplot2); library(forcats); library(patchwork)
})

# ---- Config ----
stations_15 <- c("BFL","BAB","FAT","HJO","MAE","MYV","MCE","OVE",
                 "RBL","RDD","SAC","SMF","SCK","VIS","PTV")

# Clean p01m: "T" (trace) -> 0; "null"/"" -> NA; numeric strings -> numeric
clean_p01m <- function(x) {
  x_chr <- as.character(x)
  out <- suppressWarnings(as.numeric(x_chr))
  out[is.na(out) & tolower(x_chr) %in% c("null","")] <- NA_real_
  out[toupper(x_chr) == "T"] <- 0
  out
}

# Add time fields + numeric precip; filter to 15 stations
prep_time_precip <- function(df) {
  df %>%
    filter(station %in% stations_15) %>%
    mutate(
      year  = year(local_time),
      month_num = month(local_time),
      month = factor(month.abb[month_num], levels = month.abb, ordered = TRUE),
      p01m_num = clean_p01m(p01m)
    )
}

# ---------------- Data ----------------
dust <- dust_events_cv7 %>% prep_time_precip()

# use same year window as dust so panels are comparable
yr_min <- min(dust$year, na.rm = TRUE)
yr_max <- max(dust$year, na.rm = TRUE)

climo <- asos_data19 %>%
  filter(station %in% stations_15,
         between(year(local_time), yr_min, yr_max)) %>%
  prep_time_precip()

# ---------------- Yearly totals (station-sum, then average across stations) ----
sum_yearly_by_station <- function(df, label) {
  df %>%
    group_by(station, year) %>%
    summarise(total_mm = sum(p01m_num, na.rm = TRUE), .groups = "drop") %>%
    group_by(year) %>%
    summarise(total_mm = mean(total_mm, na.rm = TRUE), .groups = "drop") %>%
    mutate(dataset = label)
}

yr_dust  <- sum_yearly_by_station(dust,  "Dust events")
yr_climo <- sum_yearly_by_station(climo, "Climatology")
yearly_totals <- bind_rows(yr_climo, yr_dust)

# ---------------- Monthly totals (station-sum per year, then average across stations & years) ----
sum_monthly_station_year <- function(df, label) {
  df %>%
    group_by(station, year, month, month_num) %>%
    summarise(total_mm = sum(p01m_num, na.rm = TRUE), .groups = "drop") %>%
    group_by(month, month_num) %>%
    summarise(total_mm = mean(total_mm, na.rm = TRUE), .groups = "drop") %>%
    mutate(dataset = label)
}

mon_dust  <- sum_monthly_station_year(dust,  "Dust events")
mon_climo <- sum_monthly_station_year(climo, "Climatology")
monthly_totals <- bind_rows(mon_climo, mon_dust) %>% arrange(month_num)

# ---------------- Plot styling ----------------
pal <- c(`Climatology` = "grey80", `Dust events` = "#2C7FB8")

base_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---------------- Plots ----------------
p_year <- ggplot(yearly_totals, aes(x = year, y = total_mm, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_fill_manual(values = pal) +
  labs(y = "Yearly total precipitation (mm)\n(station-sum, averaged across 15 stations)") +
  scale_x_continuous(breaks = seq(yr_min, yr_max, by = 2)) +
  base_theme

p_mon <- ggplot(monthly_totals, aes(x = month, y = total_mm, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_fill_manual(values = pal) +
  labs(y = "Monthly total precipitation (mm)\n(station-sum per yr, then averaged across stations & years)",
       x = NULL) +
  base_theme

p_precip <- p_year / p_mon
p_precip
# ggsave("precip_totals_dust_vs_climo.png", p_precip, width = 14, height = 8, dpi = 300)



devtools::install_github("SantanderMetGroup/convertR")
