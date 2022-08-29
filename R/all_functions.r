#' @title Scale Zero To One
#' @description Scale a numeric vector such that the minimum is 0 and maximum is 1.
#' @param .vec Numeric vector without NA's.
#' @return Numeric vector of same len of input.
scale_zero_to_one <- function(.vec){
  return((.vec - min(.vec)) / max(.vec - min(.vec) ))
}


#' @title Load and Preprocess
#' @description Loads a '.wav' file and applies several steps to identify
#' ticks.
#' @importFrom magrittr %>% %<>% not
#' @import dplyr
#' @param .path_to_wav
#' @param .lower_freqlim
#' @param .upper_freqlim
#' @return
#' @export
load_preprocess <- function(.path_to_wav,
                            .lower_freqlim = 8*10^4,
                            .upper_freqlim = 2*10^5,
                            .outlier_sensitivity = 25){

  #load file
  wave = tuneR::readWave(.path_to_wav) #To read the first file
  # plot(wave)

  #get samples per second
  samprate = wave@samp.rate

  #get total amount of samples
  sample_amt <- length(wave)

  #determine length of recording in seconds
  time_amt <- sample_amt/samprate

  #create a evenly spaced time-index vector
  timesteps <- seq(0, time_amt, by=1/samprate)
  #make sure there are as many steps as entries in wave
  timesteps <- timesteps[1:length(wave)]
  #convert to period class (slow)
  timesteps <- timesteps %>% lubridate::seconds()

  #extract amplitude information of wave, keep original info for later
  ori_wave <- wave@left

  #create copy, here, transformation will be applie
  wave <- ori_wave

  #isolate ultra-sound range in frequency-domain
  #this bypasses some of the typical recording-artifacts
  wave <- seewave::bwfilter(wave, samprate, from = .lower_freqlim, to = .upper_freqlim, output = "Wave") %>%
    .@left #get the numeric values from the wave object, again.
  # plot(wave, type='l')

  #substract running mean (trend) to minimize noise/artifacts
  #running mean is applied every 10 timesteps
  wave <-  wave - (forecast::ma(wave, 10)) %>% as.numeric()
  wave <- ifelse(wave < 0, 0, wave) #focus only on upper channel for now
  wave[is.na(wave)] <- 0 # in case there are NA's, fill
  # plot(wave, type= 'l')

  #combine info in tibble
  df <- tibble::tibble(timestamp = timesteps,
               wave = wave %>% scale_zero_to_one(),
               ori_wave = ori_wave)

  #apply inter-quartal-range outlier detection and apply score
  out <- outliers::scores(df$wave, type = 'iqr')

  #create vector (0 or 1) which determines if outlier is bigger then score 10 or not
  df$dots <- as.integer(out > .outlier_sensitivity)
  df$outlier_score <- out
  df <- df %>% mutate(
    pos = ifelse(dots==1, dots + wave, NA), # create a value to identify the marked ticks in plots
  )
  df$samprate <- samprate
  return(df)
}#endfunc



#' @title Identify Click Groups
#' @description Add variable to input data frame which isolates tick groups.
#' @import tidyr
#' @import purrr
#' @param .preprocessed_df
#' @param .radius
#' @param .threshold
#' @return Same df as input df, but with grouping var.
#' @export
identify_outlier_groups <- function(.preprocessed_df,
                                    .radius = 500,#in indx
                                    .threshold = 1/(2*.radius),
                                    .padding_length = 10^-3 * 0.05, #in seconds
                                    .additional_idx_padding = 100){ #amt_idx for correction

  amt_indices <- length(.preprocessed_df$dots)

  cli::cli_progress_bar(total = amt_indices,
                        format = "Identifying Groups | {cli::pb_bar} {cli::pb_percent}")

  grouping <- c()
  dot <- 1
  for(idx in seq(1, amt_indices) ){
    left_idx <- max(idx-.radius,  0)
    right_idx <- min(idx+.radius, amt_indices)
    if(mean(.preprocessed_df$dots[left_idx:right_idx]) >= .threshold){
      grouping[idx] <- dot
      dot <- dot+1
    }else{
      grouping[idx] <- 0
    }
    cli::cli_progress_update()
  }#enfor
  .preprocessed_df$grouping <- grouping

  #######################################

  cli::cli_progress_bar(total = amt_indices,
                        format = "Marking Groups | {cli::pb_bar} {cli::pb_percent}")

  tmp <- c(rep(0, amt_indices))
  group <- 1
  for(idx in seq(1, amt_indices)){
    right_idx <- ifelse(idx+1 >= length(grouping), length(grouping), idx+1)
    if(grouping[right_idx] > grouping[idx]){
      tmp[idx] <- group
    }else{
      group <- group + 1
    }
    cli::cli_progress_update()
  }#endfor
  tmp <- tmp/max(tmp, na.rm = T)

  # plot(tmp)

  .preprocessed_df$grouping <- tmp

  #######################################
  # for each group, get the the click's peak
  inds <- .preprocessed_df %>%
    # mutate(index = row_number()) %>%
    select(grouping, ori_wave, timestamp) %>%
    group_by(grouping) %>%
    nest(data = c(timestamp, ori_wave)) %>%
    ungroup() %>%
    mutate(peak = map_int(data, ~max(.x$ori_wave)),
           peak_ts = map2_dbl(peak, data, ~.y$timestamp[.y$ori_wave==.x] %>% .[1])
    ) %>%
    filter(peak > 0, grouping > 0, !is.na(grouping)) %>%
    pull(peak_ts)

  #locate the indices of the local peaks
  inds <- which(.preprocessed_df$timestamp %in% inds)

  .preprocessed_df <- .preprocessed_df %>% mutate(
    click_pos = NA
  )

  # match the position of the local peaks of the transformed wave to the original wave
  .preprocessed_df$click_pos[inds] <- .preprocessed_df$ori_wave[inds]
  .preprocessed_df <- .preprocessed_df %>%
    select(timestamp, ori_wave, click_pos, samprate)
  return(.preprocessed_df)
}#endfunc

calc_power_and_spectra <- function(.grouping_df,
                                   .padding_len = 0.05){
  #how many timesteps correspond to 1 ms ?
  rec_len <- max(.grouping_df$timestamp) %>% ceiling()  #len in seconds of recording
  amt_indices <- nrow(.grouping_df)
  padding_len <- 10^-3 * .padding_len #desired padding in s
  time_amt_per_idx <- rec_len/amt_indices
  padding_idx_amt <-  padding_len/time_amt_per_idx
  padding_idx_amt <- as.integer(padding_idx_amt) #force int for a nice index

  clicks_overall_timeidx_r <- seq(1, padding_idx_amt) *time_amt_per_idx #convert idx to time
  clicks_overall_timeidx_l <- rev(clicks_overall_timeidx_r) *-1 #mirror to the left
  clicks_overall_timeidx <- c(clicks_overall_timeidx_l, 0, clicks_overall_timeidx_r) * 1000 #assemble time-axis

  .grouping_df$cenered_wave <- .grouping_df$ori_wave %>% scale() %>% .[,1]
  samprate <- .grouping_df$samprate[1]

  calcs_df <- tibble(
    click_peak_ts = .grouping_df$timestamp[!is.na(.grouping_df$click_pos)],
    lower_idx = which(!is.na(.grouping_df$click_pos)) - padding_idx_amt,
    upper_idx = which(!is.na(.grouping_df$click_pos)) + padding_idx_amt
  ) %>% mutate(
    clicks = map2(lower_idx, upper_idx, ~tibble(click_timestamp = clicks_overall_timeidx,
                                                click_signals = .grouping_df$ori_wave[.x:.y],
                                                click_signals_centered = .grouping_df$cenered_wave[.x:.y])),

    # clicks = map2(lower_idx, upper_idx, ~df$ori_wave[.x:.y]),
    click_peak_ts = as.numeric(click_peak_ts),
    pow_spec_dens = map(clicks, ~seewave::spec(.x$click_signals, f = samprate, PSD=TRUE, dB="max0", plot = F) ),
    pow_spec_dens_frec = map(pow_spec_dens, ~.x[,1] ),
    pow_spec_dens_power = map(pow_spec_dens, ~.x[,2])
  ) %>%
    select(-lower_idx, -upper_idx, -pow_spec_dens)

  return(calcs_df)
}

gen_frec_spec_plot <- function(.calcs_df){
  .calcs_df %>% select(-clicks) %>%
    unnest(c(pow_spec_dens_frec, pow_spec_dens_power)) %>%
    ggplot(aes(x=pow_spec_dens_frec, y=pow_spec_dens_power, color=click_peak_ts, group=click_peak_ts)) +
    labs(title = glue::glue("Frequency spectrum of the clicks from file {fs::path_ext_remove(basename(firstfile))}"),
         subtitle = glue::glue("Recorded the {fs::file_info(firstfile)$modification_time} | Recording length: {rec_len} seconds"),
         caption = glue::glue("Sampling rate: {samprate}. Y-axis: 'max0' dB.")) +
    xlab('Frequency [Hz]') + ylab('Amplitude [dBA]') +
    geom_line(size=1, alpha=0.6) +
    scale_colour_viridis_c(option = "plasma", name = "Recording moment:") +
    hrbrthemes::theme_ft_rc()  +
    theme(plot.background = element_rect(fill = "gray3"),
          panel.background = element_rect(fill = "gray3", color = 'gray3')
    ) +
    theme(legend.position = "bottom")
}


gen_pow_clicks_plots <- function(.calcs_df){
  .calcs_df %>% select(-pow_spec_dens_frec, -pow_spec_dens_power) %>%
    unnest(clicks) %>%
    ggplot(aes(x=click_timestamp, y=click_signals_centered, color=click_peak_ts, group=click_peak_ts)) +
    labs(title = glue::glue("Clicks from file '{fs::path_ext_remove(basename(firstfile))}'"),
         subtitle = glue::glue("Recorded the {fs::file_info(firstfile)$modification_time} | Recording length: {rec_len} seconds"),
         caption = glue::glue("Click peaks windowed by {padding_len} seconds. Click power has been z-scored to enable overplotting.")) +
    xlab('Time [ms]') + ylab('Power [z-scored dBA]') +
    geom_line(size=1, alpha=0.6) +
    scale_colour_viridis_c(option = "plasma", name = "Recording moment:") +
    hrbrthemes::theme_ft_rc()  +
    geom_hline(yintercept = 0, color = 'red', size=0.5, alpha=1, linetype=2) +
    geom_vline(xintercept = 0, color = 'red', size=0.5, alpha=1, linetype=2) +
    scale_x_continuous(labels = scales::comma) +
    theme(plot.background = element_rect(fill = "gray3"),
          panel.background = element_rect(fill = "gray3", color = 'gray3')
    ) +
    theme(legend.position = "bottom")
}


interval_analysis <- function(.grouping_df){
  tibble(peak_times = .grouping_df$timestamp[!is.na(.grouping_df$click_pos)]%>% as.numeric())  %>%
    mutate(
      click_intervals = c(0, (diff(peak_times) %>% as.numeric())),
      interval_changes = c(0, (diff(click_intervals) %>% as.numeric())),
      click = seq(1, nrow(.))
    )
}

gen_interval_analysis_plot <- function(.interval_df){
  interval_analysis_df %>%
    pivot_longer(-click) %>%
    ggplot(aes(x=click, y=value)) +
    # labs(title = glue::glue("Clicks from file '{fs::path_ext_remove(basename(firstfile))}'"),
    #      subtitle = glue::glue("Recorded the {fs::file_info(firstfile)$modification_time} | Recording length: {rec_len} seconds"),
    #      caption = glue::glue("Click peaks windowed by {padding_len} seconds. Click power has been z-scored to enable overplotting.")) +
    # xlab('Time [ms]') + ylab('Power [z-scored dBA]') +
    geom_line(size=1, alpha=0.6, color='white') +
    hrbrthemes::theme_ft_rc()  +
    theme(plot.background = element_rect(fill = "gray3"),
          panel.background = element_rect(fill = "gray3", color = 'gray3')
    ) + facet_grid(vars(name), scales = 'free_y') +
    geom_hline(yintercept = 0, color = 'red', size=0.5, alpha=1, linetype=2)
}
