# Now load up the metadata
meta <- read.table('tables/original_sample_metadata.csv', sep = '\t', header = T, quote = '', comment.char = '#')

# Turn the day of week into days since feeding
meta$days.since.feed <- numeric(nrow(meta))
meta$days.since.feed[meta$day.of.week == 'Monday'] <- 3
meta$days.since.feed[meta$day.of.week == 'Tuesday'] <- 1
meta$days.since.feed[meta$day.of.week == 'Wednesday'] <- 2
meta$days.since.feed[meta$day.of.week == 'Thursday'] <- 1
meta$days.since.feed[meta$day.of.week == 'Friday'] <- 2
meta$days.since.feed[meta$day.of.week == 'Saturday'] <- 1
meta$days.since.feed[meta$day.of.week == 'Sunday'] <- 2

new.meta <- meta[, c('sample', 'date', 'age', 'sex', 'box.uses', 'days.since.feed')]
colnames(new.meta)[colnames(new.meta) == 'date'] <- 'collection.date'

write.table(new.meta, file = 'tables/sample_metadata.csv', sep = '\t', row.names = F, quote = F)

