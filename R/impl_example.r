
wd = getwd()
files = list.files(path = fs::path(wd, 'wav_feb19'), pattern = ".wav", full.names = T) #To list the files
firstfile <- files[23]  #23 is very simple #6 for debugging

df <-  clickpeakfindr::load_preprocess(.path_to_wav = firstfile)

df <- clickpeakfindr::identify_outlier_groups(df)

outliersplot <- gen_outliers_plot(df)
outliersplot

calcs_df <- calc_power_and_spectra(df)

frecspecplot <- gen_frec_spec_plot(calcs_df)
frecspecplot

peakpowplot <- gen_pow_clicks_plots(calcs_df)
peakpowplot

interval_analysis_df <- interval_analysis(df)
interval_plot <- gen_interval_analysis_plot(interval_analysis_df)
interval_plot
