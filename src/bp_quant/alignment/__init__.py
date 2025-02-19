STRINGENT_BOWTIE_PARAMS = ["--no-mixed",
                           "--dpad", "0",
                           "--gbar", "99999999",
                           "--mp", "1,1",
                           "--np", "1",
                           "--score-min", "L,0,-0.01"]

STRINGENT_STAR_PARAMS = ["-outFilterMismatchNoverReadLmax", "0.3",
                         "--scoreDelOpen", "-2",
                         "--scoreInsOpen", "-2",
                         "--scoreDelBase", "-2",
                         "--scoreInsBase", "-2"]