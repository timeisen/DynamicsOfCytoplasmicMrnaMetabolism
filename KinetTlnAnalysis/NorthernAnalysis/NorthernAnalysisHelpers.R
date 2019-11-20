library(tidyverse)
library(scales)

DataProcess <- function(tailsHT, tailsIntensity){

	##These lines remove the negative values only, but don't further process.
	# tailsHT <- read_tsv("figures/other/otherV9/NM_007907.txt") %>%
	# 	group_by(sample) %>%
	# 	mutate(count = ifelse(count < 0, 0, count)) %>% 
	# 	# filter(tail_length < 250) %>% #exclude 250
	# 	mutate(intensity = count / sum(count)) %>%
	# 	ungroup()


	#Additional tibble
	additionalHTdata = tibble(
		tail_length = rep(251:400, each = length(unique(tailsHT$sample))), 
			sample = rep(unique(tailsHT$sample), 150),
			intensity = rep(0, 150 * length(unique(tailsHT$sample))))

	##These lines average the 250mer aross all values between 250 and 400
	tailsHT <- tailsHT %>%
		bind_rows(additionalHTdata) %>%
		group_by(sample) %>%
		# mutate(count = ifelse(count < 0, 0, count)) %>% 
		mutate(count250 = ifelse(tail_length == 250, count/150, 0)) %>%
		mutate(count = ifelse(tail_length >= 250, sum(count250, na.rm = TRUE), count)) %>%
		# filter(tail_length < 250) %>% #exclude 250
		mutate(intensity = count / sum(count, na.rm = TRUE)) %>%
		ungroup()
	
	tailsIntensity <- tailsIntensity %>%
		filter(!is.na(tail_length) & tail_length > 0) %>%
		# add_row(tailsIntensity, tail_length = rep(250, length(unique(tailsIntensity$sample))), 
			# sample = unique(tailsIntensity$sample), 
			# intensity = rep(0, length(unique(tailsIntensity$sample)))) %>%
		group_by(sample) %>%
		# mutate(longtailCount = ifelse(tail_length >= 250, sum(intensity, na.rm = TRUE), 0)) %>%
		# mutate(intensity = ifelse(tail_length == 250, longtailCount,intensity)) %>%
		# filter(tail_length < 250) %>% #Excluding the 250, after all that work. 
		mutate(intensity_norm = intensity / sum(intensity,na.rm = TRUE)) %>%
		ungroup() 
	
	sLevels = rev(c("SS","8 h","4 h","2 h","1 h","40 min"))
	tailsIntensity$sample <- factor(tailsIntensity$sample, levels = sLevels)
	tailsHT$sample <- factor(tailsHT$sample, levels = sLevels)
	
	## binning
	b <- c(-Inf, 1:79 * 5, Inf) #vector of breaks
	bnames <- c(0, 1:79 * 5) #one fewer than b
	
	tailsIntensity <- mutate(tailsIntensity, tail_length.cat = cut(tailsIntensity$tail_length, breaks = b, labels = bnames)) %>% 
		group_by(tail_length.cat, sample) %>% 
		mutate(intensity_norm_bin = mean(intensity_norm)) %>%
		ungroup() %>%
		mutate(tail_length = ifelse(tail_length < 20 & sample != 'SS', NA, tail_length))
	
	tailsHT <- mutate(tailsHT, tail_length.cat = cut(tailsHT$tail_length, breaks = b, labels = bnames)) %>% 
		group_by(tail_length.cat, sample) %>% 
		mutate(intensity_norm_bin = mean(intensity)) %>%
		ungroup() %>%
		mutate(tail_length = ifelse(tail_length < 20 & sample != 'SS', NA, tail_length))
	
	hline = tibble(
		y = 0,
		yend = 0, 
		x = c(rep(20,5),0),
		xend = 400,
		sample = unique(tailsHT$sample))
	
	p13 <- ggplot(tailsIntensity,aes(x = tail_length, y = intensity_norm_bin, color = sample)) +
	geom_step(size=0.5,direction="hv") +
	scale_colour_manual(values=cols)  +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand=c(0,0),limits=c(0,400), breaks = c(0,20,100,200,300,400)) +
	theme_tim() +
	geom_step(data = tailsHT, aes(x = tail_length, y = intensity_norm_bin, color = sample), linetype = '11') + 
	facet_wrap(~sample, ncol = 1, scales = 'free_x') + 
	geom_segment(data = hline, aes(x = x, y = y, xend = xend, yend = yend), color  = 'black', size = 0.2) + 
	theme(
	legend.position="none",
	axis.line.x=element_blank())

	p12 <- ggplot(tailsIntensity,aes(x = tail_length, y = intensity_norm_bin, color = sample)) +
	geom_step(size=0.5,direction="hv") +
	scale_colour_manual(values=cols)  +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand=c(0,0),limits=c(0,400), breaks = c(0,20,100,200,300,350)) +
	theme_tim() +
	geom_step(data = tailsHT, aes(x = tail_length, y = intensity_norm_bin, color = sample), linetype = '11') + 
	facet_wrap(~sample, ncol = 1, scales = 'free_x') + 
	geom_segment(data = hline, aes(x = x, y = y, xend = xend, yend = yend), color  = 'black', size = 0.2) + 
	theme(
	legend.position="none",
	axis.line.x=element_blank())

	ggsave(plot = p12, file = 'figures/other/otherV9/NorthernCompPostdocTalk.pdf', width = 2, height = 6)

	return(p13)
}

DataProcessThbs1 <- function(tailsHT, tailsIntensity){

	##These lines remove the negative values only, but don't further process.
	# tailsHT <- read_tsv("figures/other/otherV9/NM_007907.txt") %>%
	# 	group_by(sample) %>%
	# 	mutate(count = ifelse(count < 0, 0, count)) %>% 
	# 	# filter(tail_length < 250) %>% #exclude 250
	# 	mutate(intensity = count / sum(count)) %>%
	# 	ungroup()


	#Additional tibble
	additionalHTdata = tibble(
		tail_length = rep(251:350, each = length(unique(tailsIntensity$sample))), 
			sample = rep(unique(tailsIntensity$sample), 100),
			intensity = rep(0, 100 * length(unique(tailsIntensity$sample))))

	##These lines average the 250mer aross all values between 250 and 400
	tailsHT <- tailsHT %>%
		bind_rows(additionalHTdata) %>%
		group_by(sample) %>%
		mutate(count = ifelse(count < 0, 0, count)) %>% 
		mutate(count250 = ifelse(tail_length == 250, count/100, 0)) %>%
		mutate(count = ifelse(tail_length >= 250, sum(count250, na.rm = TRUE), count)) %>%
		# filter(tail_length < 250) %>% #exclude 250
		mutate(intensity = count / sum(count, na.rm = TRUE)) %>%
		ungroup()
	tailsHT %>% filter(tail_length == 250)  %>% print()
	
	tailsIntensity <- tailsIntensity %>%
		filter(!is.na(tail_length) & tail_length > 0) %>%
		# add_row(tailsIntensity, tail_length = rep(250, length(unique(tailsIntensity$sample))), 
			# sample = unique(tailsIntensity$sample), 
			# intensity = rep(0, length(unique(tailsIntensity$sample)))) %>%
		group_by(sample) %>%
		# mutate(longtailCount = ifelse(tail_length >= 250, sum(intensity, na.rm = TRUE), 0)) %>%
		# mutate(intensity = ifelse(tail_length == 250, longtailCount,intensity)) %>%
		# filter(tail_length < 250) %>% #Excluding the 250, after all that work. 
		mutate(intensity_norm = intensity / sum(intensity,na.rm = TRUE)) %>%
		ungroup() 
	tailsIntensity %>% distinct(sample, .keep_all = TRUE) %>% print()
	
	sLevels = rev(c("SS","8 h","4 h","2 h","1 h","40 min"))
	tailsIntensity$sample <- factor(tailsIntensity$sample, levels = sLevels)
	tailsHT$sample <- factor(tailsHT$sample, levels = sLevels)
	
	## binning
	b <- c(-Inf, 1:79 * 5, Inf) #vector of breaks
	bnames <- c(0, 1:79 * 5) #one fewer than b
	
	tailsIntensity <- mutate(tailsIntensity, tail_length.cat = cut(tailsIntensity$tail_length, breaks = b, labels = bnames)) %>% 
		group_by(tail_length.cat, sample) %>% 
		mutate(intensity_norm_bin = mean(intensity_norm)) %>%
		ungroup() %>%
		mutate(tail_length = ifelse(tail_length < 20 & sample != 'SS', NA, tail_length))

	tailsHT <- mutate(tailsHT, tail_length.cat = cut(tailsHT$tail_length, breaks = b, labels = bnames)) %>% 
		group_by(tail_length.cat, sample) %>% 
		mutate(intensity_norm_bin = mean(intensity)) %>%
		ungroup() %>%
		mutate(tail_length = ifelse(tail_length < 20 & sample != 'SS', NA, tail_length))
	
	hline = tibble(
		y = 0,
		yend = 0, 
		x = c(rep(20,5),0),
		xend = 350,
		sample = unique(tailsHT$sample))
	
	p13 <- ggplot(tailsIntensity,aes(x = tail_length, y = intensity_norm_bin, color = sample)) +
	geom_step(size=0.5,direction="hv") +
	scale_colour_manual(values=cols)  +
	scale_y_continuous(expand = c(0,0)) +
	scale_x_continuous(expand=c(0,0),limits=c(0,400), breaks = c(0,20,100,200,300,350)) +
	theme_tim() +
	geom_step(data = tailsHT, aes(x = tail_length, y = intensity_norm_bin, color = sample), linetype = '11') + 
	facet_wrap(~sample, ncol = 1, scales = 'free') + 
	geom_segment(data = hline, aes(x = x, y = y, xend = xend, yend = yend), color  = 'black', size = 0.2) + 
	theme(
	legend.position="none",
	axis.line.x=element_blank())



	return(p13)
}