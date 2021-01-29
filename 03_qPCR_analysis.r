library(plotrix)

hk.genes <- c('GMOY003952', 'GMOY010976')
age.genes <- c('GMOY000749', 'GMOY001603', 'GMOY002920', 'GMOY003090', 'GMOY003371', 'GMOY003588', 'GMOY005053', 'GMOY005321', 'GMOY009908', 'GMOY011979')
genes <- c('GMOY000749', 'GMOY001603', 'GMOY002920', 'GMOY003090', 'GMOY003371', 'GMOY003588', 'GMOY003952', 'GMOY005053', 'GMOY005321', 'GMOY009908', 'GMOY010976', 'GMOY011979')

# Set the info for each plate
plate.info <- list(
	plate1 = list(Path = 'qPCR/190731/eric-tsetse-p1-190731 - Tabular Results.txt', Samples = c('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'NTC'), Replicate = 1, Genes = genes),
	plate2 = list(Path = 'qPCR/190731/eric-tsetse-p2-190731 - Tabular Results.txt', Samples = c('T8', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'NTC'), Replicate = 1, Genes = genes),
	plate3 = list(Path = 'qPCR/190731/eric-tsetse-p3-190731 - Tabular Results.txt', Samples = c('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'NTC'), Replicate = 2, Genes = genes),
	plate4 = list(Path = 'qPCR/190731/eric-tsetse-p4-190731 - Tabular Results.txt', Samples = c('T8', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'NTC'), Replicate = 2, Genes = genes),
	plate5 = list(Path = 'qPCR/190801/eric-tsetse-p5-190801 - Tabular Results.txt', Samples = c('T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'NTC'), Replicate = 1, Genes = genes),
	plate6 = list(Path = 'qPCR/190801/eric-tsetse-p6-190801 - Tabular Results.txt', Samples = c('T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34', 'NTC'), Replicate = 1, Genes = genes),
	plate7 = list(Path = 'qPCR/190801/eric-tsetse-p7-190801 - Tabular Results.txt', Samples = c('T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'NTC'), Replicate = 2, Genes = genes),
	plate8 = list(Path = 'qPCR/190801/eric-tsetse-p8-190801 - Tabular Results.txt', Samples = c('T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34', 'NTC'), Replicate = 2, Genes = genes),
	plate9 = list(Path = 'qPCR/190802/eric-tsetse-p9-190802 - Tabular Results.txt', Samples = c('T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'NTC'), Replicate = 1, Genes = genes),
	plate10 = list(Path = 'qPCR/190802/eric-tsetse-p10-190802 - Tabular Results.txt', Samples = c('T42', 'T43', 'T44', 'T45', 'T46', 'T47', 'T48', 'NTC'), Replicate = 1, Genes = genes),
	plate11 = list(Path = 'qPCR/190802/eric-tsetse-p11-190802 - Tabular Results.txt', Samples = c('T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'NTC'), Replicate = 2, Genes = genes),
	plate12 = list(Path = 'qPCR/190802/eric-tsetse-p12-190802 - Tabular Results.txt', Samples = c('T42', 'T43', 'T44', 'T45', 'T46', 'T47', 'T48', 'NTC'), Replicate = 2, Genes = genes),
	plate13 = list(Path = 'qPCR/190807/eric-tsetse-p13-190807 - Tabular Results.txt', Samples = c('T49', 'T50', 'T51', 'T52', 'T53', 'T54', 'T101', 'NTC'), Replicate = 1, Genes = genes),
	plate14 = list(Path = 'qPCR/190807/eric-tsetse-p14-190807 - Tabular Results.txt', Samples = c('T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'NTC'), Replicate = 1, Genes = genes),
	# There was a mistake in plates 15 and 16, where two columns (ie: primers) were swapped, so deal with this here
	plate15 = list(Path = 'qPCR/190807/eric-tsetse-p15-190807 - Tabular Results.txt', Samples = c('T49', 'T50', 'T51', 'T52', 'T53', 'T54', 'T101', 'NTC'), Replicate =2, Genes = genes[c(1:7, 9, 8, 10:12)]),
	plate16 = list(Path = 'qPCR/190807/eric-tsetse-p16-190807 - Tabular Results.txt', Samples = c('T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'NTC'), Replicate = 2, Genes = genes[c(1:7, 9, 8, 10:12)]),
	plate17 = list(Path = 'qPCR/190808/eric-tsetse-p17-190808 - Tabular Results.txt', Samples = c('T109', 'T110', 'T111', 'T112', 'T113', 'T114', 'T115', 'NTC'), Replicate = 1, Genes = genes),
	plate18 = list(Path = 'qPCR/190808/eric-tsetse-p18-190808 - Tabular Results.txt', Samples = c('T116', 'T117', 'T118', 'T119', 'T120', 'T121', 'T122', 'NTC'), Replicate = 1, Genes = genes),
	plate19 = list(Path = 'qPCR/190808/eric-tsetse-p19-190808 - Tabular Results.txt', Samples = c('T109', 'T110', 'T111', 'T112', 'T113', 'T114', 'T115', 'NTC'), Replicate = 2, Genes = genes),
	plate20 = list(Path = 'qPCR/190808/eric-tsetse-p20-190808 - Tabular Results.txt', Samples = c('T116', 'T117', 'T118', 'T119', 'T120', 'T121', 'T122', 'NTC'), Replicate = 2, Genes = genes),
	# There was a mistake in plates 21 and 22, where two columns (ie: primers) were swapped, so deal with this here
	plate21 = list(Path = 'qPCR/190809/eric-tsetse-p21-190809 - Tabular Results.txt', Samples = c('T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T129', 'NTC'), Replicate = 1, Genes = genes[c(1:7, 9, 8, 10:12)]),
	plate22 = list(Path = 'qPCR/190809/eric-tsetse-p22-190809 - Tabular Results.txt', Samples = c('T130', 'T131', 'T132', 'T133', 'T134', 'T135', 'T136', 'NTC'), Replicate = 1, Genes = genes[c(1:7, 9, 8, 10:12)]),
	plate23 = list(Path = 'qPCR/190809/eric-tsetse-p23-190809 - Tabular Results.txt', Samples = c('T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T129', 'NTC'), Replicate = 2, Genes = genes),
	plate24 = list(Path = 'qPCR/190809/eric-tsetse-p24-190809 - Tabular Results.txt', Samples = c('T130', 'T131', 'T132', 'T133', 'T134', 'T135', 'T136', 'NTC'), Replicate = 2, Genes = genes),
	plate25 = list(Path = 'qPCR/190813/eric-tsetse-p25-190813 - Tabular Results.txt', Samples = c('T137', 'T138', 'T139', 'T140', 'T141', 'T142', 'T143', 'NTC'), Replicate = 1, Genes = genes),
	plate26 = list(Path = 'qPCR/190813/eric-tsetse-p26-190813 - Tabular Results.txt', Samples = c('T144', 'T145', 'T146', 'T147', 'T148', 'T149', 'T150', 'NTC'), Replicate = 1, Genes = genes),
	plate27 = list(Path = 'qPCR/190813/eric-tsetse-p27-190813 - Tabular Results.txt', Samples = c('T137', 'T138', 'T139', 'T140', 'T141', 'T142', 'T143', 'NTC'), Replicate = 2, Genes = genes),
	plate28 = list(Path = 'qPCR/190813/eric-tsetse-p28-190813 - Tabular Results.txt', Samples = c('T144', 'T145', 'T146', 'T147', 'T148', 'T149', 'T150', 'NTC'), Replicate = 2, Genes = genes),
	plate29 = list(Path = 'qPCR/190814/eric-tsetse-p29-190814 - Tabular Results.txt', Samples = c('T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'NTC'), Replicate = 1, Genes = genes),
	plate30 = list(Path = 'qPCR/190814/eric-tsetse-p30-190814 - Tabular Results.txt', Samples = c('T158', 'T159', 'T160', 'T161', 'T162', 'T163', 'T164', 'NTC'), Replicate = 1, Genes = genes),
	plate31 = list(Path = 'qPCR/190814/eric-tsetse-p31-190814 - Tabular Results.txt', Samples = c('T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'NTC'), Replicate = 2, Genes = genes),
	plate32 = list(Path = 'qPCR/190814/eric-tsetse-p32-190814 - Tabular Results.txt', Samples = c('T158', 'T159', 'T160', 'T161', 'T162', 'T163', 'T164', 'NTC'), Replicate = 2, Genes = genes),
	plate33 = list(Path = 'qPCR/190814/eric-tsetse-p33-190814 - Tabular Results.txt', Samples = c('T165', 'T166', 'T167', 'T168', 'T169', 'T170', 'T171', 'NTC'), Replicate = 1, Genes = genes),
	plate34 = list(Path = 'qPCR/190814/eric-tsetse-p34-190814 - Tabular Results.txt', Samples = c('T172', 'T173', 'T174', 'T175', 'T176', 'T177', 'T178', 'NTC'), Replicate = 1, Genes = genes),
	plate35 = list(Path = 'qPCR/190814/eric-tsetse-p35-190814 - Tabular Results.txt', Samples = c('T165', 'T166', 'T167', 'T168', 'T169', 'T170', 'T171', 'NTC'), Replicate = 2, Genes = genes),
	plate36 = list(Path = 'qPCR/190814/eric-tsetse-p36-190814 - Tabular Results.txt', Samples = c('T172', 'T173', 'T174', 'T175', 'T176', 'T177', 'T178', 'NTC'), Replicate = 2, Genes = genes),
	plate37 = list(Path = 'qPCR/190815/eric-tsetse-p37-190815 - Tabular Results.txt', Samples = c('T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T185', 'NTC'), Replicate = 1, Genes = genes),
	plate38 = list(Path = 'qPCR/190815/eric-tsetse-p38-190815 - Tabular Results.txt', Samples = c('T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T192', 'NTC'), Replicate = 1, Genes = genes),
	plate39 = list(Path = 'qPCR/190815/eric-tsetse-p39-190815 - Tabular Results.txt', Samples = c('T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T185', 'NTC'), Replicate = 2, Genes = genes),
	plate40 = list(Path = 'qPCR/190815/eric-tsetse-p40-190815 - Tabular Results.txt', Samples = c('T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T192', 'NTC'), Replicate = 2, Genes = genes),
	plate41 = list(Path = 'qPCR/190816/eric-tsetse-p41-190816 - Tabular Results.txt', Samples = c('T193', 'T194', 'T195', 'T196', 'T197', 'T198', 'T199', 'NTC'), Replicate = 1, Genes = genes),
	plate42 = list(Path = 'qPCR/190816/eric-tsetse-p42-190816 - Tabular Results.txt', Samples = c('T200', 'T201', 'T202', 'T203', 'T204', 'T205', 'T206', 'NTC'), Replicate = 1, Genes = genes),
	plate43 = list(Path = 'qPCR/190816/eric-tsetse-p43-190816 - Tabular Results.txt', Samples = c('T193', 'T194', 'T195', 'T196', 'T197', 'T198', 'T199', 'NTC'), Replicate = 2, Genes = genes),
	plate44 = list(Path = 'qPCR/190816/eric-tsetse-p44-190816 - Tabular Results.txt', Samples = c('T200', 'T201', 'T202', 'T203', 'T204', 'T205', 'T206', 'NTC'), Replicate = 2, Genes = genes),
	plate45 = list(Path = 'qPCR/190816/eric-tsetse-p45-190816 - Tabular Results.txt', Samples = c('T207', 'T208', 'T209b', 'T210a', 'T210b', 'T211', 'T212', 'NTC'), Replicate = 1, Genes = genes),
	plate46 = list(Path = 'qPCR/190816/eric-tsetse-p46-190816 - Tabular Results.txt', Samples = c('T213', 'T214', 'T215', 'T216', 'T217', 'T218', 'T219', 'NTC'), Replicate = 1, Genes = genes),
	plate47 = list(Path = 'qPCR/190816/eric-tsetse-p47-190816 - Tabular Results.txt', Samples = c('T207', 'T208', 'T209b', 'T210a', 'T210b', 'T211', 'T212', 'NTC'), Replicate = 2, Genes = genes),
	plate48 = list(Path = 'qPCR/190816/eric-tsetse-p48-190816 - Tabular Results.txt', Samples = c('T213', 'T214', 'T215', 'T216', 'T217', 'T218', 'T219', 'NTC'), Replicate = 2, Genes = genes),
	plate49 = list(Path = 'qPCR/190902/eric-tsetse-p49-190902 - Tabular Results.txt', Samples = c('T220', 'T221', 'T222', 'T223', 'T224', 'T225', 'T226', 'NTC'), Replicate = 1, Genes = genes),
	plate50 = list(Path = 'qPCR/190902/eric-tsetse-p50-190902 - Tabular Results.txt', Samples = c('T227', 'T228', 'T229', 'T230', 'T232', 'T233', 'T234', 'NTC'), Replicate = 1, Genes = genes),
	plate51 = list(Path = 'qPCR/190902/eric-tsetse-p51-190902 - Tabular Results.txt', Samples = c('T220', 'T221', 'T222', 'T223', 'T224', 'T225', 'T226', 'NTC'), Replicate = 2, Genes = genes),
	# There was a mistake on plate 52, the primers were loaded in reverse order, so deal with this here
	plate52 = list(Path = 'qPCR/190902/eric-tsetse-p52-190902 - Tabular Results.txt', Samples = c('T227', 'T228', 'T229', 'T230', 'T232', 'T233', 'T234', 'NTC'), Replicate = 2, Genes = genes[seq(length(genes), 1)]),
	plate53 = list(Path = 'qPCR/190903/eric-tsetse-p53-190903 - Tabular Results.txt', Samples = c('T235', 'T236', 'T237', 'T238', 'T239', 'T240', 'T241', 'NTC'), Replicate = 1, Genes = genes),
	plate54 = list(Path = 'qPCR/190903/eric-tsetse-p54-190903 - Tabular Results.txt', Samples = c('T242', 'T243', 'T244', 'T245', 'T246', 'T247', 'T248', 'NTC'), Replicate = 1, Genes = genes),
	plate55 = list(Path = 'qPCR/190903/eric-tsetse-p55-190903 - Tabular Results.txt', Samples = c('T235', 'T236', 'T237', 'T238', 'T239', 'T240', 'T241', 'NTC'), Replicate = 2, Genes = genes),
	plate56 = list(Path = 'qPCR/190903/eric-tsetse-p56-190903 - Tabular Results.txt', Samples = c('T242', 'T243', 'T244', 'T245', 'T246', 'T247', 'T248', 'NTC'), Replicate = 2, Genes = genes),
	plate57 = list(Path = 'qPCR/190904/eric-tsetse-p57-190904 - Tabular Results.txt', Samples = c('T249', 'T250', 'T251', 'T252', 'T253', 'T254', 'T255', 'NTC'), Replicate = 1, Genes = genes),
	plate58 = list(Path = 'qPCR/190904/eric-tsetse-p58-190904 - Tabular Results.txt', Samples = c('T256', 'T257', 'T258', 'T259', 'T260', 'T261', 'T262', 'NTC'), Replicate = 1, Genes = genes),
	plate59 = list(Path = 'qPCR/190904/eric-tsetse-p59-190904 - Tabular Results.txt', Samples = c('T249', 'T250', 'T251', 'T252', 'T253', 'T254', 'T255', 'NTC'), Replicate = 2, Genes = genes),
	# There was a mistake in plate 60, where two columns (ie: primers) were swapped, so deal with this here
	plate60 = list(Path = 'qPCR/190904/eric-tsetse-p60-190904 - Tabular Results.txt', Samples = c('T256', 'T257', 'T258', 'T259', 'T260', 'T261', 'T262', 'NTC'), Replicate = 2, Genes = genes[c(1:3, 5, 4, 6:12)]),
	plate61 = list(Path = 'qPCR/190904/eric-tsetse-p61-190904 - Tabular Results.txt', Samples = c('T263', 'T264', 'T265', 'T266', 'T267', 'T268', 'T269', 'NTC'), Replicate = 1, Genes = genes),
	plate62 = list(Path = 'qPCR/190904/eric-tsetse-p62-190904 - Tabular Results.txt', Samples = c('T270', 'T271', 'T272', 'T273', 'T274', 'T275', 'T276', 'NTC'), Replicate = 1, Genes = genes),
	plate63 = list(Path = 'qPCR/190904/eric-tsetse-p63-190904 - Tabular Results.txt', Samples = c('T263', 'T264', 'T265', 'T266', 'T267', 'T268', 'T269', 'NTC'), Replicate = 2, Genes = genes),
	# There was a mistake in plate 64, where two columns (ie: primers) were swapped, so deal with this here
	plate64 = list(Path = 'qPCR/190904/eric-tsetse-p64-190904 - Tabular Results.txt', Samples = c('T270', 'T271', 'T272', 'T273', 'T274', 'T275', 'T276', 'NTC'), Replicate = 2, Genes = genes[c(1:2, 4, 3, 5:12)]),
	plate65 = list(Path = 'qPCR/190905/eric-tsetse-p65-190905 - Tabular Results.txt', Samples = c('T277', 'T278', 'T279', 'T280', 'T281', 'T282', 'T283', 'NTC'), Replicate = 1, Genes = genes),
	plate66 = list(Path = 'qPCR/190905/eric-tsetse-p66-190905 - Tabular Results.txt', Samples = c('T284', 'T285', 'T286', 'T287', 'T288', 'T289', 'T290', 'NTC'), Replicate = 1, Genes = genes),
	plate67 = list(Path = 'qPCR/190905/eric-tsetse-p67-190905 - Tabular Results.txt', Samples = c('T277', 'T278', 'T279', 'T280', 'T281', 'T282', 'T283', 'NTC'), Replicate = 2, Genes = genes),
	plate68 = list(Path = 'qPCR/190905/eric-tsetse-p68-190905 - Tabular Results.txt', Samples = c('T284', 'T285', 'T286', 'T287', 'T288', 'T289', 'T290', 'NTC'), Replicate = 2, Genes = genes),
	plate69 = list(Path = 'qPCR/190913/eric-tsetse-p69-190913 - Tabular Results.txt', Samples = c('T291', 'T292', 'T293', 'T294', 'T295', 'T296', 'T297', 'NTC'), Replicate = 1, Genes = genes),
	plate70 = list(Path = 'qPCR/190913/eric-tsetse-p70-190913 - Tabular Results.txt', Samples = c('T298', 'T299', 'T300', 'T301', 'T302', 'T303', 'T304', 'NTC'), Replicate = 1, Genes = genes),
	plate71 = list(Path = 'qPCR/190913/eric-tsetse-p71-190913 - Tabular Results.txt', Samples = c('T291', 'T292', 'T293', 'T294', 'T295', 'T296', 'T297', 'NTC'), Replicate = 2, Genes = genes),
	plate72 = list(Path = 'qPCR/190913/eric-tsetse-p72-190913 - Tabular Results.txt', Samples = c('T298', 'T299', 'T300', 'T301', 'T302', 'T303', 'T304', 'NTC'), Replicate = 2, Genes = genes),
	plate73 = list(Path = 'qPCR/190913/eric-tsetse-p73-190913 - Tabular Results.txt', Samples = c('T305', 'T306', 'T307', 'T308', 'T309', 'T310', 'T311', 'NTC'), Replicate = 1, Genes = genes),
	plate74 = list(Path = 'qPCR/190913/eric-tsetse-p74-190913 - Tabular Results.txt', Samples = c('T312', 'T313', 'T314', 'T315', 'T316', 'T317', 'T318', 'NTC'), Replicate = 1, Genes = genes),
	plate75 = list(Path = 'qPCR/190913/eric-tsetse-p75-190913 - Tabular Results.txt', Samples = c('T305', 'T306', 'T307', 'T308', 'T309', 'T310', 'T311', 'NTC'), Replicate = 2, Genes = genes),
	plate76 = list(Path = 'qPCR/190913/eric-tsetse-p76-190913 - Tabular Results.txt', Samples = c('T312', 'T313', 'T314', 'T315', 'T316', 'T317', 'T318', 'NTC'), Replicate = 2, Genes = genes),
	plate77 = list(Path = 'qPCR/190916/eric-tsetse-p77-190916 - Tabular Results.txt', Samples = c('T319', 'T320', 'T321', 'T322', 'T323', 'T324', 'T325', 'NTC'), Replicate = 1, Genes = genes),
	plate78 = list(Path = 'qPCR/190916/eric-tsetse-p78-190916 - Tabular Results.txt', Samples = c('T326', 'T327', 'T328', 'T329', 'T330', 'T331', 'T332', 'NTC'), Replicate = 1, Genes = genes),
	plate79 = list(Path = 'qPCR/190916/eric-tsetse-p79-190916 - Tabular Results.txt', Samples = c('T319', 'T320', 'T321', 'T322', 'T323', 'T324', 'T325', 'NTC'), Replicate = 2, Genes = genes),
	plate80 = list(Path = 'qPCR/190916/eric-tsetse-p80-190916 - Tabular Results.txt', Samples = c('T326', 'T327', 'T328', 'T329', 'T330', 'T331', 'T332', 'NTC'), Replicate = 2, Genes = genes),
	plate81 = list(Path = 'qPCR/190916/eric-tsetse-p81-190916 - Tabular Results.txt', Samples = c('T333', 'T334', 'T335', 'T336', 'T337', 'T338', 'T339', 'NTC'), Replicate = 1, Genes = genes),
	plate82 = list(Path = 'qPCR/190916/eric-tsetse-p82-190916 - Tabular Results.txt', Samples = c('T340', 'T341', 'T342', 'T343', 'T344', 'T345', 'T346', 'NTC'), Replicate = 1, Genes = genes),
	plate83 = list(Path = 'qPCR/190916/eric-tsetse-p83-190916 - Tabular Results.txt', Samples = c('T333', 'T334', 'T335', 'T336', 'T337', 'T338', 'T339', 'NTC'), Replicate = 2, Genes = genes),
	plate84 = list(Path = 'qPCR/190916/eric-tsetse-p84-190916 - Tabular Results.txt', Samples = c('T340', 'T341', 'T342', 'T343', 'T344', 'T345', 'T346', 'NTC'), Replicate = 2, Genes = genes),
	plate85 = list(Path = 'qPCR/190917/eric-tsetse-p85-190917 - Tabular Results.txt', Samples = c('T347', 'T348', 'T349', 'T350', 'T351', 'T352', 'T353', 'NTC'), Replicate = 1, Genes = genes),
	plate86 = list(Path = 'qPCR/190917/eric-tsetse-p86-190917 - Tabular Results.txt', Samples = c('T354', 'T355', 'T356', 'T357', 'T358', 'T359', 'T360', 'NTC'), Replicate = 1, Genes = genes),
	plate87 = list(Path = 'qPCR/190917/eric-tsetse-p87-190917 - Tabular Results.txt', Samples = c('T347', 'T348', 'T349', 'T350', 'T351', 'T352', 'T353', 'NTC'), Replicate = 2, Genes = genes),
	plate88 = list(Path = 'qPCR/190917/eric-tsetse-p88-190917 - Tabular Results.txt', Samples = c('T354', 'T355', 'T356', 'T357', 'T358', 'T359', 'T360', 'NTC'), Replicate = 2, Genes = genes),
	plate89 = list(Path = 'qPCR/190919/eric-tsetse-p89-190919 - Tabular Results.txt', Samples = c('T361', 'T362', 'T363', 'T364', 'T365', 'T366', 'T367', 'NTC'), Replicate = 1, Genes = genes),
	plate90 = list(Path = 'qPCR/190919/eric-tsetse-p90-190919 - Tabular Results.txt', Samples = c('T368', 'T369', 'T370', 'T371', 'T372', 'T373', 'T374', 'NTC'), Replicate = 1, Genes = genes),
	plate91 = list(Path = 'qPCR/190919/eric-tsetse-p91-190919 - Tabular Results.txt', Samples = c('T361', 'T362', 'T363', 'T364', 'T365', 'T366', 'T367', 'NTC'), Replicate = 2, Genes = genes),
	plate92 = list(Path = 'qPCR/190919/eric-tsetse-p92-190919 - Tabular Results.txt', Samples = c('T368', 'T369', 'T370', 'T371', 'T372', 'T373', 'T374', 'NTC'), Replicate = 2, Genes = genes),
	plate93 = list(Path = 'qPCR/190919/eric-tsetse-p93-190919 - Tabular Results.txt', Samples = c('T375', 'T376', 'T377', 'T378', 'T379', 'T380', 'T381', 'NTC'), Replicate = 1, Genes = genes),
	plate94 = list(Path = 'qPCR/190919/eric-tsetse-p94-190919 - Tabular Results.txt', Samples = c('T382', 'T383', 'T384', 'T385', 'T386', 'T387', 'T388', 'NTC'), Replicate = 1, Genes = genes),
	plate95 = list(Path = 'qPCR/190919/eric-tsetse-p95-190919 - Tabular Results.txt', Samples = c('T375', 'T376', 'T377', 'T378', 'T379', 'T380', 'T381', 'NTC'), Replicate = 2, Genes = genes),
	plate96 = list(Path = 'qPCR/190919/eric-tsetse-p96-190919 - Tabular Results.txt', Samples = c('T382', 'T383', 'T384', 'T385', 'T386', 'T387', 'T388', 'NTC'), Replicate = 2, Genes = genes),
	plate97 = list(Path = 'qPCR/190920/eric-tsetse-p97-190920 - Tabular Results.txt', Samples = c('T389', 'T390', 'T391', 'T392', 'T393', 'T394', 'T395', 'NTC'), Replicate = 1, Genes = genes),
	plate98 = list(Path = 'qPCR/190920/eric-tsetse-p98-190920 - Tabular Results.txt', Samples = c('T396', 'T397', 'T398', 'T399', 'T400', 'T401', 'T402', 'NTC'), Replicate = 1, Genes = genes),
	plate99 = list(Path = 'qPCR/190920/eric-tsetse-p99-190920 - Tabular Results.txt', Samples = c('T389', 'T390', 'T391', 'T392', 'T393', 'T394', 'T395', 'NTC'), Replicate = 2, Genes = genes),
	# There was a mistake in plate 100, where two columns (ie: primers) were swapped, so deal with this here
	plate100 = list(Path = 'qPCR/190920/eric-tsetse-p100-190920 - Tabular Results.txt', Samples = c('T396', 'T397', 'T398', 'T399', 'T400', 'T401', 'T402', 'NTC'), Replicate = 2, Genes = genes[c(1:2, 4, 3, 5:12)]),
	plate101 = list(Path = 'qPCR/190920/eric-tsetse-p101-190920 - Tabular Results.txt', Samples = c('T403', 'T404', 'T405', 'T406', 'T407', 'T408', 'T409', 'NTC'), Replicate = 1, Genes = genes),
	plate102 = list(Path = 'qPCR/190920/eric-tsetse-p102-190920 - Tabular Results.txt', Samples = c('T410', 'T411', 'T412', 'T413', 'T414', 'T415', 'T416', 'NTC'), Replicate = 1, Genes = genes),
	plate103 = list(Path = 'qPCR/190920/eric-tsetse-p103-190920 - Tabular Results.txt', Samples = c('T403', 'T404', 'T405', 'T406', 'T407', 'T408', 'T409', 'NTC'), Replicate = 2, Genes = genes),
	plate104 = list(Path = 'qPCR/190920/eric-tsetse-p104-190920 - Tabular Results.txt', Samples = c('T410', 'T411', 'T412', 'T413', 'T414', 'T415', 'T416', 'NTC'), Replicate = 2, Genes = genes),
	plate105 = list(Path = 'qPCR/190923/eric-tsetse-p105-190923 - Tabular Results.txt', Samples = c('T417', 'T418', 'T419', 'T420', 'T421', 'T422', 'T423', 'NTC'), Replicate = 1, Genes = genes),
	plate106 = list(Path = 'qPCR/190923/eric-tsetse-p106-190923 - Tabular Results.txt', Samples = c('T424', 'T425', 'T426', 'T427', 'T428', 'T429', 'T430', 'NTC'), Replicate = 1, Genes = genes),
	plate107 = list(Path = 'qPCR/190923/eric-tsetse-p107-190923 - Tabular Results.txt', Samples = c('T417', 'T418', 'T419', 'T420', 'T421', 'T422', 'T423', 'NTC'), Replicate = 2, Genes = genes),
	plate108 = list(Path = 'qPCR/190923/eric-tsetse-p108-190923 - Tabular Results.txt', Samples = c('T424', 'T425', 'T426', 'T427', 'T428', 'T429', 'T430', 'NTC'), Replicate = 2, Genes = genes),
	plate109 = list(Path = 'qPCR/190923/eric-tsetse-p109-190923 - Tabular Results.txt', Samples = c('T431', 'T432', 'T433', 'T434', 'T435', 'T436', 'T437', 'NTC'), Replicate = 1, Genes = genes),
	plate110 = list(Path = 'qPCR/190923/eric-tsetse-p110-190923 - Tabular Results.txt', Samples = c('T438', 'T439', 'T440', 'T441', 'T442', 'T443', 'T444', 'NTC'), Replicate = 1, Genes = genes),
	plate111 = list(Path = 'qPCR/190923/eric-tsetse-p111-190923 - Tabular Results.txt', Samples = c('T431', 'T432', 'T433', 'T434', 'T435', 'T436', 'T437', 'NTC'), Replicate = 2, Genes = genes),
	plate112 = list(Path = 'qPCR/190923/eric-tsetse-p112-190923 - Tabular Results.txt', Samples = c('T438', 'T439', 'T440', 'T441', 'T442', 'T443', 'T444', 'NTC'), Replicate = 2, Genes = genes),
	plate113 = list(Path = 'qPCR/190924/eric-tsetse-p113-190924 - Tabular Results.txt', Samples = c('T445', 'T446', 'T447', 'T448', 'T449', 'T450', 'T451', 'NTC'), Replicate = 1, Genes = genes),
	plate114 = list(Path = 'qPCR/190924/eric-tsetse-p114-190924 - Tabular Results.txt', Samples = c('T452', 'T453', 'T454', 'T455', 'T456', 'T457', 'T458', 'NTC'), Replicate = 1, Genes = genes),
	plate115 = list(Path = 'qPCR/190924/eric-tsetse-p115-190924 - Tabular Results.txt', Samples = c('T445', 'T446', 'T447', 'T448', 'T449', 'T450', 'T451', 'NTC'), Replicate = 2, Genes = genes),
	plate116 = list(Path = 'qPCR/190924/eric-tsetse-p116-190924 - Tabular Results.txt', Samples = c('T452', 'T453', 'T454', 'T455', 'T456', 'T457', 'T458', 'NTC'), Replicate = 2, Genes = genes),
	plate117 = list(Path = 'qPCR/190925/eric-tsetse-p117-190925 - Tabular Results.txt', Samples = c('T459', 'T460', 'T461', 'T462', 'T463', 'T464', 'T465', 'NTC'), Replicate = 1, Genes = genes),
	plate118 = list(Path = 'qPCR/190925/eric-tsetse-p118-190925 - Tabular Results.txt', Samples = c('T466', 'T467', 'T468', 'T469', 'T470', 'T471', 'T472', 'NTC'), Replicate = 1, Genes = genes),
	plate119 = list(Path = 'qPCR/190925/eric-tsetse-p119-190925 - Tabular Results.txt', Samples = c('T459', 'T460', 'T461', 'T462', 'T463', 'T464', 'T465', 'NTC'), Replicate = 2, Genes = genes),
	plate120 = list(Path = 'qPCR/190925/eric-tsetse-p120-190925 - Tabular Results.txt', Samples = c('T466', 'T467', 'T468', 'T469', 'T470', 'T471', 'T472', 'NTC'), Replicate = 2, Genes = genes),
	# For some reason, column 11 went wrong on either plates 121&122 or plates 123&124, so the two technical 
	# replicates gave different results. Since column 11 is a housekeeping gene, we have to rerun the whole plates. 
	plate121 = list(Path = 'qPCR/190926/eric-tsetse-p121-190926 - Tabular Results.txt', Samples = c('T473', 'T474', 'T475', 'T476', 'T477', 'T478', 'T479', 'NTC'), Replicate = 1, Genes = genes),
	plate122 = list(Path = 'qPCR/190926/eric-tsetse-p122-190926 - Tabular Results.txt', Samples = c('T480', 'T481', 'T482', 'T483', 'T484', 'T485', 'T486', 'NTC'), Replicate = 1, Genes = genes),
	plate123 = list(Path = 'qPCR/190926/eric-tsetse-p123-190926 - Tabular Results.txt', Samples = c('T473', 'T474', 'T475', 'T476', 'T477', 'T478', 'T479', 'NTC'), Replicate = 2, Genes = genes),
	plate124 = list(Path = 'qPCR/190926/eric-tsetse-p124-190926 - Tabular Results.txt', Samples = c('T480', 'T481', 'T482', 'T483', 'T484', 'T485', 'T486', 'NTC'), Replicate = 2, Genes = genes),
	plate125 = list(Path = 'qPCR/190927/eric-tsetse-p125-190927 - Tabular Results.txt', Samples = c('T487', 'T488', 'T489', 'T490', 'T491', 'T492', 'T493', 'NTC'), Replicate = 1, Genes = genes),
	plate126 = list(Path = 'qPCR/190927/eric-tsetse-p126-190927 - Tabular Results.txt', Samples = c('T494', 'T495', 'T496', 'T497', 'T498', 'T499', 'T500', 'NTC'), Replicate = 1, Genes = genes),
	plate127 = list(Path = 'qPCR/190927/eric-tsetse-p127-190927 - Tabular Results.txt', Samples = c('T487', 'T488', 'T489', 'T490', 'T491', 'T492', 'T493', 'NTC'), Replicate = 2, Genes = genes),
	plate128 = list(Path = 'qPCR/190927/eric-tsetse-p128-190927 - Tabular Results.txt', Samples = c('T494', 'T495', 'T496', 'T497', 'T498', 'T499', 'T500', 'NTC'), Replicate = 2, Genes = genes),
	plate129 = list(Path = 'qPCR/190927/eric-tsetse-p129-190927 - Tabular Results.txt', Samples = c('T501', 'T502', 'T503', 'T504', 'T505', 'T506', 'T507', 'NTC'), Replicate = 1, Genes = genes),
	plate130 = list(Path = 'qPCR/190927/eric-tsetse-p130-190927 - Tabular Results.txt', Samples = c('T508', 'T509', 'T510', 'T511', 'T512', 'T513', 'T514', 'NTC'), Replicate = 1, Genes = genes),
	plate131 = list(Path = 'qPCR/190927/eric-tsetse-p131-190927 - Tabular Results.txt', Samples = c('T501', 'T502', 'T503', 'T504', 'T505', 'T506', 'T507', 'NTC'), Replicate = 2, Genes = genes),
	plate132 = list(Path = 'qPCR/190927/eric-tsetse-p132-190927 - Tabular Results.txt', Samples = c('T508', 'T509', 'T510', 'T511', 'T512', 'T513', 'T514', 'NTC'), Replicate = 2, Genes = genes),
	plate133 = list(Path = 'qPCR/190930/eric-tsetse-p133-190930 - Tabular Results.txt', Samples = c('T515', 'T516', 'T517', 'T518', 'T519', 'T520', 'T521', 'NTC'), Replicate = 1, Genes = genes),
	plate134 = list(Path = 'qPCR/190930/eric-tsetse-p134-190930 - Tabular Results.txt', Samples = c('T522', 'T523', 'T524', 'T525', 'T526', 'T527', 'T528', 'NTC'), Replicate = 1, Genes = genes),
	plate135 = list(Path = 'qPCR/190930/eric-tsetse-p135-190930 - Tabular Results.txt', Samples = c('T515', 'T516', 'T517', 'T518', 'T519', 'T520', 'T521', 'NTC'), Replicate = 2, Genes = genes),
	plate136 = list(Path = 'qPCR/190930/eric-tsetse-p136-190930 - Tabular Results.txt', Samples = c('T522', 'T523', 'T524', 'T525', 'T526', 'T527', 'T528', 'NTC'), Replicate = 2, Genes = genes),
	plate137 = list(Path = 'qPCR/190930/eric-tsetse-p137-190930 - Tabular Results.txt', Samples = c('T529', 'T530', 'T531', 'T532', 'T533', 'T534', 'T535', 'NTC'), Replicate = 1, Genes = genes),
	plate138 = list(Path = 'qPCR/190930/eric-tsetse-p138-190930 - Tabular Results.txt', Samples = c('T536', 'T537', 'T538', 'T539', 'T540', 'T541', 'T542', 'NTC'), Replicate = 1, Genes = genes),
	plate139 = list(Path = 'qPCR/190930/eric-tsetse-p139-190930 - Tabular Results.txt', Samples = c('T529', 'T530', 'T531', 'T532', 'T533', 'T534', 'T535', 'NTC'), Replicate = 2, Genes = genes),
	plate140 = list(Path = 'qPCR/190930/eric-tsetse-p140-190930 - Tabular Results.txt', Samples = c('T536', 'T537', 'T538', 'T539', 'T540', 'T541', 'T542', 'NTC'), Replicate = 2, Genes = genes),
	plate141 = list(Path = 'qPCR/191001/eric-tsetse-p141-191001 - Tabular Results.txt', Samples = c('T543', 'T544', 'T545', 'T546', 'T547', 'T548', 'T549', 'NTC'), Replicate = 1, Genes = genes),
	plate142 = list(Path = 'qPCR/191001/eric-tsetse-p142-191001 - Tabular Results.txt', Samples = c('T550', 'T231', 'NTC'), Replicate = 1, Genes = genes),
	plate143 = list(Path = 'qPCR/191001/eric-tsetse-p143-191001 - Tabular Results.txt', Samples = c('T543', 'T544', 'T545', 'T546', 'T547', 'T548', 'T549', 'NTC'), Replicate = 2, Genes = genes),
	plate144 = list(Path = 'qPCR/191001/eric-tsetse-p144-191001 - Tabular Results.txt', Samples = c('T550', 'T231', 'NTC'), Replicate = 2, Genes = genes),
	plate145 = list(Path = 'qPCR/191003/eric-tsetse-p145-191001 - Tabular Results.txt', Samples = c('T473', 'T474', 'T475', 'T476', 'T477', 'T478', 'T479', 'NTC'), Replicate = 3, Genes = genes),
	plate146 = list(Path = 'qPCR/191003/eric-tsetse-p146-191001 - Tabular Results.txt', Samples = c('T480', 'T481', 'T482', 'T483', 'T484', 'T485', 'T486', 'NTC'), Replicate = 3, Genes = genes)
)

# Write a function to load the data from the plates
load.plate <- function(plate.name){
	info.list <- plate.info[[plate.name]]
	# Do some checks
	if (any(names(info.list) != c('Path', 'Samples', 'Replicate', 'Genes')))
		stop('"info.list" does not contain the right entry names.')
	if (length(info.list[['Samples']]) != 8)
		cat('Warning: "Samples" entry for plate', plate.name, 'does not contain 8 values.\n')
	if (length(info.list[['Genes']]) != 12)
		stop('"Genes" entry should contain exactly 12 values.')
	this.plate <- read.table(info.list[['Path']], sep = '\t', header = T, row.names = 1, na.strings = 'No Cq', quote = '', strip.white = T, col.names = c('Well', 'Ct', 'Tm', 'Gene'))
	# If there are fewer than 8 samples, we only keep as many rows as there are samples
	s <- info.list[['Samples']]
	if (length(s) < 8)
		this.plate <- this.plate[1:(12*length(s)), ]
	this.plate$Gene <- info.list[['Genes']]
	this.plate$Sample <- factor(rep(s, each = 12))
	this.plate$Replicate <- info.list[['Replicate']]
	this.plate
}

plates <- plates.combined <- plates.normalised <- plates.raw <- list()
for (p in names(plate.info))
	plates.raw[[p]] <- load.plate(p)

check.hk <- function(plate){
	bob <- subset(plate, Gene %in% hk.genes & Sample != 'NTC')
	if (any(is.na(bob$Ct))) return(5)
	else return(0)
}

for (p in names(plates.raw)){
	if (check.hk(plates.raw[[p]]) == 5)
		print(p)
}

# Write a function to normalise the Ct values within a plate by the housekeeping genes
normalise.plate <- function(plate){
	output.plate <- plate
	# Replace missing values of age genes with 40
	output.plate$Ct[is.na(output.plate$Ct) & (output.plate$Gene %in% age.genes)] <- 40
	output.plate$Normalised.Ct <- numeric(nrow(output.plate))
	for (this.sample in levels(output.plate$Sample)){
		these.rows <- output.plate$Sample == this.sample
		output.plate$Normalised.Ct[these.rows] <- mean(output.plate$Ct[these.rows & output.plate$Gene %in% hk.genes]) - output.plate$Ct[these.rows]
	}
	output.plate
}

for (p in names(plates.raw))
	plates.normalised[[p]] <- normalise.plate(plates.raw[[p]])


# Write a function that combines two tables by taking the mean of the Ct, Tm and normalised Ct
combine.tables <- function(list.of.tables, tablename1, tablename2, with.plot = F, plot.diff.thresh = 1){
	table1 <- list.of.tables[[tablename1]]
	table2 <- list.of.tables[[tablename2]]
	if (any(dim(table1) != dim(table2)) | any(rownames(table1) != rownames(table2)) | any(colnames(table1) != colnames(table2)))
		stop('The two tables need to have identical rows and columns.')
	table1$Well <- rownames(table1)
	table2$Well <- rownames(table2)
	rownames(table1) <- interaction(table1$Gene, table1$Sample)
	rownames(table2) <- interaction(table2$Gene, table2$Sample)
	output.table <- data.frame(row.names = rownames(table1))
	# If the genes in table1 and table2 are not in the same order, we re-arrange table2 to match table1, 
	# giving a warning
	if (any(table1$Gene != table2$Gene)){
		cat('Tables ', tablename1, ' and ', tablename2, ' do not have genes in the same order. This makes me sad. Having to do things the hard way.\n', sep = '')
		# Fix the values in table2 so that they match the order in table 1
		table2$Gene <- table1$Gene
		table2$Ct <- table2[rownames(table1), 'Ct']
		table2$Tm <- table2[rownames(table1), 'Tm']
		table2$Normalised.Ct <- table2[rownames(table1), 'Normalised.Ct']
	}
	# Record the values of interest. These are the Ct, the Tm, and the difference in Ct between the two replicates.
	# Since some of the samples needed to be rerun, already prepare empty columns for Ct3 and Norm.Ct3 (otherwise
	# we won't be able to rbind the samples that were rerun with the ones that weren't).
	output.table$Ct1 <- table1$Ct
	output.table$Ct2 <- table2$Ct
	output.table$Ct3 <- NA
	output.table$Ct <- (table1$Ct + table2$Ct) / 2
	output.table$Tm <- (table1$Tm + table2$Tm) / 2
	output.table$Norm.Ct1 <- table1$Normalised.Ct
	output.table$Norm.Ct2 <- table2$Normalised.Ct
	output.table$Norm.Ct3 <- NA
	output.table$Normalised.Ct <- (table1$Normalised.Ct + table2$Normalised.Ct) / 2
	output.table$Norm.Ct.diff <- table1$Normalised.Ct - table2$Normalised.Ct
	output.table[, c('Gene', 'Sample', 'Replicate')] <- table1[, c('Gene', 'Sample', 'Replicate')]
	if (with.plot){
		x11()
		plot(table1$Normalised.Ct, table2$Normalised.Ct, xlab = tablename1, ylab = tablename2, type = 'n')
		Ct.diff.too.big <- abs(output.table$Norm.Ct.diff) >= plot.diff.thresh
		text(table1$Normalised.Ct, table2$Normalised.Ct, labels = table1$Well, col = c('blue', 'red')[Ct.diff.too.big + 1])
		abline(0,1)
	}
	output.table
}

plates.combined[['plates1_3']] <- combine.tables(plates.normalised, 'plate1', 'plate3')
plates.combined[['plates2_4']] <- combine.tables(plates.normalised, 'plate2', 'plate4')
plates.combined[['plates5_7']] <- combine.tables(plates.normalised, 'plate5', 'plate7')
plates.combined[['plates6_8']] <- combine.tables(plates.normalised, 'plate6', 'plate8')
plates.combined[['plates9_11']] <- combine.tables(plates.normalised, 'plate9', 'plate11')
plates.combined[['plates10_12']] <- combine.tables(plates.normalised, 'plate10', 'plate12')
plates.combined[['plates13_15']] <- combine.tables(plates.normalised, 'plate13', 'plate15')
plates.combined[['plates14_16']] <- combine.tables(plates.normalised, 'plate14', 'plate16')
plates.combined[['plates17_19']] <- combine.tables(plates.normalised, 'plate17', 'plate19')
plates.combined[['plates18_20']] <- combine.tables(plates.normalised, 'plate18', 'plate20')
plates.combined[['plates21_23']] <- combine.tables(plates.normalised, 'plate21', 'plate23')
plates.combined[['plates22_24']] <- combine.tables(plates.normalised, 'plate22', 'plate24')
plates.combined[['plates25_27']] <- combine.tables(plates.normalised, 'plate25', 'plate27')
plates.combined[['plates26_28']] <- combine.tables(plates.normalised, 'plate26', 'plate28')
plates.combined[['plates29_31']] <- combine.tables(plates.normalised, 'plate29', 'plate31')
plates.combined[['plates30_32']] <- combine.tables(plates.normalised, 'plate30', 'plate32')
plates.combined[['plates33_35']] <- combine.tables(plates.normalised, 'plate33', 'plate35')
plates.combined[['plates34_36']] <- combine.tables(plates.normalised, 'plate34', 'plate36')
plates.combined[['plates37_39']] <- combine.tables(plates.normalised, 'plate37', 'plate39')
plates.combined[['plates38_40']] <- combine.tables(plates.normalised, 'plate38', 'plate40')
plates.combined[['plates41_43']] <- combine.tables(plates.normalised, 'plate41', 'plate43')
plates.combined[['plates42_44']] <- combine.tables(plates.normalised, 'plate42', 'plate44')
plates.combined[['plates45_47']] <- combine.tables(plates.normalised, 'plate45', 'plate47')
plates.combined[['plates46_48']] <- combine.tables(plates.normalised, 'plate46', 'plate48')
plates.combined[['plates49_51']] <- combine.tables(plates.normalised, 'plate49', 'plate51')
plates.combined[['plates50_52']] <- combine.tables(plates.normalised, 'plate50', 'plate52')
plates.combined[['plates53_55']] <- combine.tables(plates.normalised, 'plate53', 'plate55')
plates.combined[['plates54_56']] <- combine.tables(plates.normalised, 'plate54', 'plate56')
plates.combined[['plates57_59']] <- combine.tables(plates.normalised, 'plate57', 'plate59')
plates.combined[['plates58_60']] <- combine.tables(plates.normalised, 'plate58', 'plate60')
plates.combined[['plates61_63']] <- combine.tables(plates.normalised, 'plate61', 'plate63')
plates.combined[['plates62_64']] <- combine.tables(plates.normalised, 'plate62', 'plate64')
plates.combined[['plates65_67']] <- combine.tables(plates.normalised, 'plate65', 'plate67')
plates.combined[['plates66_68']] <- combine.tables(plates.normalised, 'plate66', 'plate68')
plates.combined[['plates69_71']] <- combine.tables(plates.normalised, 'plate69', 'plate71')
plates.combined[['plates70_72']] <- combine.tables(plates.normalised, 'plate70', 'plate72')
plates.combined[['plates73_75']] <- combine.tables(plates.normalised, 'plate73', 'plate75')
plates.combined[['plates74_76']] <- combine.tables(plates.normalised, 'plate74', 'plate76')
plates.combined[['plates77_79']] <- combine.tables(plates.normalised, 'plate77', 'plate79')
plates.combined[['plates78_80']] <- combine.tables(plates.normalised, 'plate78', 'plate80')
plates.combined[['plates81_83']] <- combine.tables(plates.normalised, 'plate81', 'plate83')
plates.combined[['plates82_84']] <- combine.tables(plates.normalised, 'plate82', 'plate84')
plates.combined[['plates85_87']] <- combine.tables(plates.normalised, 'plate85', 'plate87')
plates.combined[['plates86_88']] <- combine.tables(plates.normalised, 'plate86', 'plate88')
plates.combined[['plates89_91']] <- combine.tables(plates.normalised, 'plate89', 'plate91')
plates.combined[['plates90_92']] <- combine.tables(plates.normalised, 'plate90', 'plate92')
plates.combined[['plates93_95']] <- combine.tables(plates.normalised, 'plate93', 'plate95')
plates.combined[['plates94_96']] <- combine.tables(plates.normalised, 'plate94', 'plate96')
plates.combined[['plates97_99']] <- combine.tables(plates.normalised, 'plate97', 'plate99')
plates.combined[['plates98_100']] <- combine.tables(plates.normalised, 'plate98', 'plate100')
plates.combined[['plates101_103']] <- combine.tables(plates.normalised, 'plate101', 'plate103')
plates.combined[['plates102_104']] <- combine.tables(plates.normalised, 'plate102', 'plate104')
plates.combined[['plates105_107']] <- combine.tables(plates.normalised, 'plate105', 'plate107')
plates.combined[['plates106_108']] <- combine.tables(plates.normalised, 'plate106', 'plate108')
plates.combined[['plates109_111']] <- combine.tables(plates.normalised, 'plate109', 'plate111')
plates.combined[['plates110_112']] <- combine.tables(plates.normalised, 'plate110', 'plate112')
plates.combined[['plates113_115']] <- combine.tables(plates.normalised, 'plate113', 'plate115')
plates.combined[['plates114_116']] <- combine.tables(plates.normalised, 'plate114', 'plate116')
plates.combined[['plates117_119']] <- combine.tables(plates.normalised, 'plate117', 'plate119')
plates.combined[['plates118_120']] <- combine.tables(plates.normalised, 'plate118', 'plate120')
plates.combined[['plates121_123']] <- combine.tables(plates.normalised, 'plate121', 'plate123')
plates.combined[['plates122_124']] <- combine.tables(plates.normalised, 'plate122', 'plate124')
plates.combined[['plates125_127']] <- combine.tables(plates.normalised, 'plate125', 'plate127')
plates.combined[['plates126_128']] <- combine.tables(plates.normalised, 'plate126', 'plate128')
plates.combined[['plates129_131']] <- combine.tables(plates.normalised, 'plate129', 'plate131')
plates.combined[['plates130_132']] <- combine.tables(plates.normalised, 'plate130', 'plate132')
plates.combined[['plates133_135']] <- combine.tables(plates.normalised, 'plate133', 'plate135')
plates.combined[['plates134_136']] <- combine.tables(plates.normalised, 'plate134', 'plate136')
plates.combined[['plates137_139']] <- combine.tables(plates.normalised, 'plate137', 'plate139')
plates.combined[['plates138_140']] <- combine.tables(plates.normalised, 'plate138', 'plate140')
plates.combined[['plates141_143']] <- combine.tables(plates.normalised, 'plate141', 'plate143')
plates.combined[['plates142_144']] <- combine.tables(plates.normalised, 'plate142', 'plate144')

# Look for samples that might need repeating. 
all.plates.combined <- do.call(rbind, plates.combined)
possible.reruns <- subset(all.plates.combined, abs(Norm.Ct.diff) > 1)

# Overall, we need to rerun:
# On plates 1-4: genes 000749 and 001603 + hk (so four columns 16 rows) = 2/3 plate
# On plates 5-8: gene 003371 + hk (so three columns 16 rows) = 1/2 plate
# On plates 121-124: whole plate (so 12 columns 16 rows) = 2 plates
# On plates 37-40: gene 000749 + hk (so three columns 16 rows) = 1/2 plate

# There was a lot of variability for gene 005321, but this was in the low expression samples and was 
# probably just noise inherent in samples with low expression.

# We create a reduced table that excludes all the variation described in the above text, leaving the 
# individual samples that failed. 
rows.to.remove <- grepl('plates1_', rownames(possible.reruns)) |
                  grepl('plates2_', rownames(possible.reruns)) |
                  grepl('plates5_', rownames(possible.reruns)) |
                  grepl('plates6_', rownames(possible.reruns)) |
                  grepl('plates121_', rownames(possible.reruns)) |
                  grepl('plates122_', rownames(possible.reruns)) |
                  grepl('plates37_', rownames(possible.reruns)) |
                  grepl('plates38_', rownames(possible.reruns))

reduced.reruns <- subset(possible.reruns[!rows.to.remove, ], Gene != 'GMOY005321')
reruns.321 <- subset(possible.reruns[!rows.to.remove, ], Gene == 'GMOY005321')

# Let's rerun all of those. 
# On gene 001603, we need to rerun the following samples:
#T47 T102 T205 T239 T369 T418 T534: That's 7 samples + NTC on three columns

# On gene 002920, we need to rerun the following samples:
#T162 T236 T382

# On gene 002908, we need to rerun the following samples:
#T162 T432 T502

# On gene 000749, we need to rerun the following sample:
#T135

# On gene 003588, we need to rerun the following sample:
#T456

thresh.lines <- function(){
	abline(0,1)
	abline(1,1, col = 'grey50')
	abline(-1,1, col = 'grey50')
}

# Right, now lets look at the samples that we repeated. 
plate.142.repeat.genes <- c('GMOY003371', hk.genes)
plate.142 <- read.table(plate.info[['plate142']][['Path']], sep = '\t', header = T, row.names = 1, na.strings = 'No Cq', quote = '', strip.white = T, col.names = c('Well', 'Ct', 'Tm', 'Gene'))
# Keep only the wells that were repeats
repeat.wells <- c(paste(rep(c('D', 'E', 'F', 'G'), each = 9), 1:9, sep = ''), paste('H', 1:6, sep = ''))
plate.142.repeats <- plate.142[repeat.wells, ]
plate.142.repeats$Gene <- plate.142.repeat.genes
s <- c('T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34')
plate.142.repeats$Sample <- factor(rep(s, each = 3))
plate.142.repeats <- normalise.plate(plate.142.repeats)
rownames(plate.142.repeats) <- interaction(plate.142.repeats$Gene, plate.142.repeats$Sample)
# Now add these values to the other plates with the same samples
plates.combined[['plates5_7']]$Ct3 <- plate.142.repeats[rownames(plates.combined[['plates5_7']]), 'Ct']
plates.combined[['plates5_7']]$Norm.Ct3 <- plate.142.repeats[rownames(plates.combined[['plates5_7']]), 'Normalised.Ct']
# Now plot the three Normalised Cts. We first plot the second (red) and third (blue) replicates against the 
# first, then we plot the first (green) and third(blue) replicates againse the second. We can see that the
# third replicate closely matches the first, but not the second replicate, so we consider that the second
# replicate was wrong and we will keep the first and third replicates, averaging the two. This applies to 
# the samples from plates 5_7 and 6_8
plates5_7 <- plates.combined[['plates5_7']]
#
x11()
par(mfrow = c(2,2))
not.na <- !is.na(plates5_7$Ct3)
these.pch <- rep(1, nrow(plates5_7[not.na, ]))
these.pch[plates5_7$Gene[not.na] %in% hk.genes] <- 3
#
plot(plates5_7$Norm.Ct1[not.na], plates5_7$Norm.Ct2[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'red', pch = these.pch, ylim = c(-1.5, 6))
points(plates5_7$Norm.Ct1[not.na], plates5_7$Norm.Ct3[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'blue', pch = these.pch)
thresh.lines()
#
plot(plates5_7$Norm.Ct2[not.na], plates5_7$Norm.Ct1[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'green', pch = these.pch)
points(plates5_7$Norm.Ct2[not.na], plates5_7$Norm.Ct3[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'blue', pch = these.pch)
thresh.lines()
# Now do plates6_8
plates.combined[['plates6_8']]$Ct3 <- plate.142.repeats[rownames(plates.combined[['plates6_8']]), 'Ct']
plates.combined[['plates6_8']]$Norm.Ct3 <- plate.142.repeats[rownames(plates.combined[['plates6_8']]), 'Normalised.Ct']
plates6_8 <- plates.combined[['plates6_8']]
#
not.na <- !is.na(plates6_8$Ct3)
these.pch <- rep(1, nrow(plates6_8[not.na, ]))
these.pch[plates6_8$Gene[not.na] %in% hk.genes] <- 3
#
plot(plates6_8$Norm.Ct1[not.na], plates6_8$Norm.Ct2[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'red', pch = these.pch, ylim = c(-1.5, 6))
points(plates6_8$Norm.Ct1[not.na], plates6_8$Norm.Ct3[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'blue', pch = these.pch)
thresh.lines()
#
plot(plates6_8$Norm.Ct2[not.na], plates6_8$Norm.Ct1[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'green', pch = these.pch)
points(plates6_8$Norm.Ct2[not.na], plates6_8$Norm.Ct3[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'blue', pch = these.pch)
thresh.lines()
# Since the first and third replicates are in agreement, we recalculate the mean Ct and normalised Ct for the gene
# that was repeated
repeated.rows.5_7 <- rownames(subset(plates.combined[['plates5_7']], Gene == plate.142.repeat.genes[1]))
plates.combined[['plates5_7']][repeated.rows.5_7, 'Ct'] <- (plates.combined[['plates5_7']][repeated.rows.5_7, 'Ct1'] + plates.combined[['plates5_7']][repeated.rows.5_7, 'Ct3']) / 2
plates.combined[['plates5_7']][repeated.rows.5_7, 'Normalised.Ct'] <- (plates.combined[['plates5_7']][repeated.rows.5_7, 'Norm.Ct1'] + plates.combined[['plates5_7']][repeated.rows.5_7, 'Norm.Ct3']) / 2
repeated.rows.6_8 <- rownames(subset(plates.combined[['plates6_8']], Gene == plate.142.repeat.genes[1]))
plates.combined[['plates6_8']][repeated.rows.6_8, 'Ct'] <- (plates.combined[['plates6_8']][repeated.rows.6_8, 'Ct1'] + plates.combined[['plates6_8']][repeated.rows.6_8, 'Ct3']) / 2
plates.combined[['plates6_8']][repeated.rows.6_8, 'Normalised.Ct'] <- (plates.combined[['plates6_8']][repeated.rows.6_8, 'Norm.Ct1'] + plates.combined[['plates6_8']][repeated.rows.6_8, 'Norm.Ct3']) / 2

# Now do the same for plates 37-40 
plate.144.repeat.genes <- c('GMOY000749', hk.genes)
plate.144 <- read.table(plate.info[['plate144']][['Path']], sep = '\t', header = T, row.names = 1, na.strings = 'No Cq', quote = '', strip.white = T, col.names = c('Well', 'Ct', 'Tm', 'Gene'))
# Keep only the wells that were repeats
repeat.wells <- paste(rep(c('D', 'E', 'F', 'G', 'H'), each = 9), 1:9, sep = '')
plate.144.repeats <- plate.144[repeat.wells, ]
plate.144.repeats$Gene <- plate.144.repeat.genes
s <- c('T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T185', 'T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T192', 'T135')
plate.144.repeats$Sample <- factor(rep(s, each = 3))
plate.144.repeats <- normalise.plate(plate.144.repeats)
rownames(plate.144.repeats) <- interaction(plate.144.repeats$Gene, plate.144.repeats$Sample)
# Now add these values to the other plates with the same samples
plates.combined[['plates37_39']]$Ct3 <- plate.144.repeats[rownames(plates.combined[['plates37_39']]), 'Ct']
plates.combined[['plates37_39']]$Norm.Ct3 <- plate.144.repeats[rownames(plates.combined[['plates37_39']]), 'Normalised.Ct']
# Now plot the three Normalised Cts. We first plot the second (red) and third (blue) replicates against the 
# first, then we plot the first (green) and third(blue) replicates againse the second. We can see that the
# third replicate closely matches the second, but not the second first, so we consider that the second
# replicate was wrong and we will keep the first and third replicates, averaging the two. This applies to 
# the samples from plates 37_39 and 38_40
plates37_39 <- plates.combined[['plates37_39']]
x11()
par(mfrow = c(2,2))
not.na <- !is.na(plates37_39$Ct3)
these.pch <- rep(1, nrow(plates37_39[not.na, ]))
these.pch[plates37_39$Gene[not.na] %in% hk.genes] <- 3
#
plot(plates37_39$Norm.Ct1[not.na], plates37_39$Norm.Ct2[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'red', pch = these.pch, xlim = c(-1.5, 6), ylim = c(-1.5, 6))
points(plates37_39$Norm.Ct1[not.na], plates37_39$Norm.Ct3[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'blue', pch = these.pch)
thresh.lines()
#
plot(plates37_39$Norm.Ct2[not.na], plates37_39$Norm.Ct1[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'green', pch = these.pch, xlim = c(-1.5, 6), ylim = c(-1.5, 6))
points(plates37_39$Norm.Ct2[not.na], plates37_39$Norm.Ct3[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'blue', pch = these.pch)
thresh.lines()
# Now do plates38_40
plates.combined[['plates38_40']]$Ct3 <- plate.144.repeats[rownames(plates.combined[['plates38_40']]), 'Ct']
plates.combined[['plates38_40']]$Norm.Ct3 <- plate.144.repeats[rownames(plates.combined[['plates38_40']]), 'Normalised.Ct']
plates38_40 <- plates.combined[['plates38_40']]
#
not.na <- !is.na(plates38_40$Ct3)
these.pch <- rep(1, nrow(plates38_40[not.na, ]))
these.pch[plates38_40$Gene[not.na] %in% hk.genes] <- 3
#
plot(plates38_40$Norm.Ct1[not.na], plates38_40$Norm.Ct2[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'red', pch = these.pch, xlim = c(-1.5, 6), ylim = c(-1.5, 6))
points(plates38_40$Norm.Ct1[not.na], plates38_40$Norm.Ct3[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'blue', pch = these.pch)
thresh.lines()
#
plot(plates38_40$Norm.Ct2[not.na], plates38_40$Norm.Ct1[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'green', pch = these.pch, xlim = c(-1.5, 6), ylim = c(-1.5, 6))
points(plates38_40$Norm.Ct2[not.na], plates38_40$Norm.Ct3[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'blue', pch = these.pch)
thresh.lines()
# Since the second and third replicates are in agreement, we recalculate the mean Ct and normalised Ct for the gene
# that was repeated
repeated.rows.37_39 <- rownames(subset(plates.combined[['plates37_39']], Gene == plate.144.repeat.genes[1]))
plates.combined[['plates37_39']][repeated.rows.37_39, 'Ct'] <- (plates.combined[['plates37_39']][repeated.rows.37_39, 'Ct2'] + plates.combined[['plates37_39']][repeated.rows.37_39, 'Ct3']) / 2
plates.combined[['plates37_39']][repeated.rows.37_39, 'Normalised.Ct'] <- (plates.combined[['plates37_39']][repeated.rows.37_39, 'Norm.Ct2'] + plates.combined[['plates37_39']][repeated.rows.37_39, 'Norm.Ct3']) / 2
repeated.rows.38_40 <- rownames(subset(plates.combined[['plates38_40']], Gene == plate.144.repeat.genes[1]))
plates.combined[['plates38_40']][repeated.rows.38_40, 'Ct'] <- (plates.combined[['plates38_40']][repeated.rows.38_40, 'Ct2'] + plates.combined[['plates38_40']][repeated.rows.38_40, 'Ct3']) / 2
plates.combined[['plates38_40']][repeated.rows.38_40, 'Normalised.Ct'] <- (plates.combined[['plates38_40']][repeated.rows.38_40, 'Norm.Ct2'] + plates.combined[['plates38_40']][repeated.rows.38_40, 'Norm.Ct3']) / 2

# For plates, 121-124, it's the second replicate (plates 123 and 124) that were wrong. So we combine
# 121 with 145 and 122 with 146
plate.145 <- read.table(plate.info[['plate145']][['Path']], sep = '\t', header = T, row.names = 1, na.strings = 'No Cq', quote = '', strip.white = T, col.names = c('Well', 'Ct', 'Tm', 'Gene'))
plate.145$Gene <- plate.info[['plate145']][['Genes']]
plate.145$Sample <- factor(rep(plate.info[['plate145']][['Samples']], each = 12))
plate.145 <- normalise.plate(plate.145)
rownames(plate.145) <- interaction(plate.145$Gene, plate.145$Sample)
plates.combined[['plates121_123']]$Ct3 <- plate.145[rownames(plates.combined[['plates121_123']]), 'Ct']
plates.combined[['plates121_123']]$Norm.Ct3 <- plate.145[rownames(plates.combined[['plates121_123']]), 'Normalised.Ct']
#
these.pch <- rep(1, nrow(plates.combined[['plates121_123']]))
these.pch[plates.combined[['plates121_123']]$Gene %in% hk.genes] <- 3
x11()
par(mfrow = c(2,2))
plot(plates.combined[['plates121_123']]$Norm.Ct1, plates.combined[['plates121_123']]$Norm.Ct2, , xlab = 'Rep1', ylab = 'Rep2/3', col = 'red', pch = these.pch, ylim = c(-1.5, 6))
points(plates.combined[['plates121_123']]$Norm.Ct1, plates.combined[['plates121_123']]$Norm.Ct3, xlab = 'Rep1', ylab = 'Rep2/3', col = 'blue', pch = these.pch)
thresh.lines()
plot(plates.combined[['plates121_123']]$Norm.Ct2, plates.combined[['plates121_123']]$Norm.Ct1, , xlab = 'Rep2', ylab = 'Rep1/3', col = 'green', pch = these.pch, ylim = c(-1.5, 6))
points(plates.combined[['plates121_123']]$Norm.Ct2, plates.combined[['plates121_123']]$Norm.Ct3, xlab = 'Rep2', ylab = 'Rep1/3', col = 'blue', pch = these.pch)
thresh.lines()
#
plate.146 <- read.table(plate.info[['plate146']][['Path']], sep = '\t', header = T, row.names = 1, na.strings = 'No Cq', quote = '', strip.white = T, col.names = c('Well', 'Ct', 'Tm', 'Gene'))
plate.146$Gene <- plate.info[['plate146']][['Genes']]
plate.146$Sample <- factor(rep(plate.info[['plate146']][['Samples']], each = 12))
plate.146 <- normalise.plate(plate.146)
rownames(plate.146) <- interaction(plate.146$Gene, plate.146$Sample)
plates.combined[['plates122_124']]$Ct3 <- plate.146[rownames(plates.combined[['plates122_124']]), 'Ct']
plates.combined[['plates122_124']]$Norm.Ct3 <- plate.146[rownames(plates.combined[['plates122_124']]), 'Normalised.Ct']
#
plot(plates.combined[['plates122_124']]$Norm.Ct1, plates.combined[['plates122_124']]$Norm.Ct2, , xlab = 'Rep1', ylab = 'Rep2/3', col = 'red', pch = these.pch, ylim = c(-1.5, 6))
points(plates.combined[['plates122_124']]$Norm.Ct1, plates.combined[['plates122_124']]$Norm.Ct3, xlab = 'Rep1', ylab = 'Rep2/3', col = 'blue', pch = these.pch)
thresh.lines()
plot(plates.combined[['plates122_124']]$Norm.Ct2, plates.combined[['plates122_124']]$Norm.Ct1, , xlab = 'Rep2', ylab = 'Rep1/3', col = 'green', pch = these.pch, ylim = c(-1.5, 6))
points(plates.combined[['plates122_124']]$Norm.Ct2, plates.combined[['plates122_124']]$Norm.Ct3, xlab = 'Rep2', ylab = 'Rep1/3', col = 'blue', pch = these.pch)
thresh.lines()
# This shows that the first and third replicates agree with each other, with the second replicate being an outlier. 
# We therefore re-calculate the Ct means using the first and third replicates. 
plates.combined[['plates121_123']]$Ct <- (plates.combined[['plates121_123']]$Ct1 + plates.combined[['plates121_123']]$Ct3)/2
plates.combined[['plates121_123']]$Normalised.Ct <- (plates.combined[['plates121_123']]$Norm.Ct1 + plates.combined[['plates121_123']]$Norm.Ct3)/2
plates.combined[['plates122_124']]$Ct <- (plates.combined[['plates122_124']]$Ct1 + plates.combined[['plates122_124']]$Ct3)/2
plates.combined[['plates122_124']]$Normalised.Ct <- (plates.combined[['plates122_124']]$Norm.Ct1 + plates.combined[['plates122_124']]$Norm.Ct3)/2


# Now plate 147, which repeated samples on plates 1-4 for genes 000703 and 001603, and some other samples
# for 001603
plate.147 <- read.table('qPCR/191004/eric-tsetse-p147-191004 - Tabular Results.txt', sep = '\t', header = T, row.names = 1, na.strings = 'No Cq', quote = '', strip.white = T, col.names = c('Well', 'Ct', 'Tm', 'Gene'))
# Group 1 contains the repeats for plates 1-4. Group 2 contains other repeats. 
group1.wells <- paste(rep(LETTERS[1:8], each = 8), 1:8, sep = '')
group2.wells <- paste(rep(LETTERS[1:7], each = 3), 9:11, sep = '')
plate.147.group1 <- plate.147[group1.wells, ]
plate.147.group2 <- plate.147[group2.wells, ]
# Let's start with the repeats of plates 1-4 first
group1.genes <- c('GMOY000749', 'GMOY001603', hk.genes)
plate.147.group1$Gene <- group1.genes
group1.s <- c('T1', 'T8', 'T2', 'T14', 'T3', 'T15', 'T4', 'T16', 'T5', 'T17', 'T6', 'T19', 'T7', 'T20', 'NTC1', 'NTC2')
plate.147.group1$Sample <- factor(rep(group1.s, each = 4))
plate.147.group1 <- normalise.plate(plate.147.group1)
rownames(plate.147.group1) <- interaction(plate.147.group1$Gene, plate.147.group1$Sample)
# Now add these values to the other plates with the same samples
plates.combined[['plates1_3']]$Ct3 <- plate.147.group1[rownames(plates.combined[['plates1_3']]), 'Ct']
plates.combined[['plates1_3']]$Norm.Ct3 <- plate.147.group1[rownames(plates.combined[['plates1_3']]), 'Normalised.Ct']
# Now plot the three Normalised Cts. We first plot the second (red) and third (blue) replicates against the 
# first, then we plot the first (green) and third(blue) replicates againse the second. 
plates1_3 <- plates.combined[['plates1_3']]
x11()
par(mfrow = c(2,2))
not.na <- !is.na(plates1_3$Ct3)
these.pch <- rep(1, nrow(plates1_3[not.na, ]))
these.pch[plates1_3$Gene[not.na] %in% hk.genes] <- 3
#
plot(plates1_3$Norm.Ct1[not.na], plates1_3$Norm.Ct2[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'red', pch = these.pch, xlim = c(-1.5, 6), ylim = c(-1.5, 6))
points(plates1_3$Norm.Ct1[not.na], plates1_3$Norm.Ct3[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'blue', pch = these.pch)
thresh.lines()
#
plot(plates1_3$Norm.Ct2[not.na], plates1_3$Norm.Ct1[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'green', pch = these.pch, xlim = c(-1.5, 6), ylim = c(-1.5, 6))
points(plates1_3$Norm.Ct2[not.na], plates1_3$Norm.Ct3[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'blue', pch = these.pch)
thresh.lines()
# Now do plates2_4
plates.combined[['plates2_4']]$Ct3 <- plate.147.group1[rownames(plates.combined[['plates2_4']]), 'Ct']
plates.combined[['plates2_4']]$Norm.Ct3 <- plate.147.group1[rownames(plates.combined[['plates2_4']]), 'Normalised.Ct']
plates2_4 <- plates.combined[['plates2_4']]
#
not.na <- !is.na(plates2_4$Ct3)
these.pch <- rep(1, nrow(plates2_4[not.na, ]))
these.pch[plates2_4$Gene[not.na] %in% hk.genes] <- 3
#
plot(plates2_4$Norm.Ct1[not.na], plates2_4$Norm.Ct2[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'red', pch = these.pch, xlim = c(-1.5, 6), ylim = c(-1.5, 6))
points(plates2_4$Norm.Ct1[not.na], plates2_4$Norm.Ct3[not.na], xlab = 'Rep1', ylab = 'Rep2/3', col = 'blue', pch = these.pch)
thresh.lines()
#
plot(plates2_4$Norm.Ct2[not.na], plates2_4$Norm.Ct1[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'green', pch = these.pch, xlim = c(-1.5, 6), ylim = c(-1.5, 6))
points(plates2_4$Norm.Ct2[not.na], plates2_4$Norm.Ct3[not.na], xlab = 'Rep2', ylab = 'Rep1/3', col = 'blue', pch = these.pch)
thresh.lines()
# Since the first and third replicates are in agreement, we recalculate the mean Ct and normalised Ct for the gene
# that was repeated
repeated.rows.1_3 <- rownames(subset(plates.combined[['plates1_3']], Gene %in% group1.genes[1:2]))
plates.combined[['plates1_3']][repeated.rows.1_3, 'Ct'] <- (plates.combined[['plates1_3']][repeated.rows.1_3, 'Ct1'] + plates.combined[['plates1_3']][repeated.rows.1_3, 'Ct3']) / 2
plates.combined[['plates1_3']][repeated.rows.1_3, 'Normalised.Ct'] <- (plates.combined[['plates1_3']][repeated.rows.1_3, 'Norm.Ct1'] + plates.combined[['plates1_3']][repeated.rows.1_3, 'Norm.Ct3']) / 2
repeated.rows.2_4 <- rownames(subset(plates.combined[['plates2_4']], Gene %in% group1.genes[1:2]))
plates.combined[['plates2_4']][repeated.rows.2_4, 'Ct'] <- (plates.combined[['plates2_4']][repeated.rows.2_4, 'Ct1'] + plates.combined[['plates2_4']][repeated.rows.2_4, 'Ct3']) / 2
plates.combined[['plates2_4']][repeated.rows.2_4, 'Normalised.Ct'] <- (plates.combined[['plates2_4']][repeated.rows.2_4, 'Norm.Ct1'] + plates.combined[['plates2_4']][repeated.rows.2_4, 'Norm.Ct3']) / 2

# Now let's look at the other repeats on this plate. 
group2.genes <- c('GMOY001603', hk.genes)
plate.147.group2$Gene <- group2.genes
group2.s <- c('T47', 'T102', 'T205', 'T239', 'T369', 'T418', 'T534')
plate.147.group2$Sample <- factor(rep(group2.s, each = 3))
plate.147.group2 <- normalise.plate(plate.147.group2)
rownames(plate.147.group2) <- interaction(plate.147.group2$Gene, plate.147.group2$Sample)
# For each repeat, let's compare with the two originals
c(plates.combined[['plates10_12']]['GMOY001603.T47', 'Norm.Ct1'], plates.combined[['plates10_12']]['GMOY001603.T47', 'Norm.Ct2'], plate.147.group2['GMOY001603.T47', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates10_12']]$Ct3 <- plate.147.group2[rownames(plates.combined[['plates10_12']]), 'Ct']
plates.combined[['plates10_12']]$Norm.Ct3 <- plate.147.group2[rownames(plates.combined[['plates10_12']]), 'Normalised.Ct']
plates.combined[['plates10_12']]['GMOY001603.T47', 'Ct'] <- (plates.combined[['plates10_12']]['GMOY001603.T47', 'Ct1'] + plates.combined[['plates10_12']]['GMOY001603.T47', 'Ct3']) / 2
plates.combined[['plates10_12']]['GMOY001603.T47', 'Normalised.Ct'] <- (plates.combined[['plates10_12']]['GMOY001603.T47', 'Norm.Ct1'] + plates.combined[['plates10_12']]['GMOY001603.T47', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates14_16']]['GMOY001603.T102', 'Norm.Ct1'], plates.combined[['plates14_16']]['GMOY001603.T102', 'Norm.Ct2'], plate.147.group2['GMOY001603.T102', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates14_16']]$Ct3 <- plate.147.group2[rownames(plates.combined[['plates14_16']]), 'Ct']
plates.combined[['plates14_16']]$Norm.Ct3 <- plate.147.group2[rownames(plates.combined[['plates14_16']]), 'Normalised.Ct']
plates.combined[['plates14_16']]['GMOY001603.T102', 'Ct'] <- (plates.combined[['plates14_16']]['GMOY001603.T102', 'Ct1'] + plates.combined[['plates14_16']]['GMOY001603.T102', 'Ct3']) / 2
plates.combined[['plates14_16']]['GMOY001603.T102', 'Normalised.Ct'] <- (plates.combined[['plates14_16']]['GMOY001603.T102', 'Norm.Ct1'] + plates.combined[['plates14_16']]['GMOY001603.T102', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates42_44']]['GMOY001603.T205', 'Norm.Ct1'], plates.combined[['plates42_44']]['GMOY001603.T205', 'Norm.Ct2'], plate.147.group2['GMOY001603.T205', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates42_44']]$Ct3 <- plate.147.group2[rownames(plates.combined[['plates42_44']]), 'Ct']
plates.combined[['plates42_44']]$Norm.Ct3 <- plate.147.group2[rownames(plates.combined[['plates42_44']]), 'Normalised.Ct']
plates.combined[['plates42_44']]['GMOY001603.T205', 'Ct'] <- (plates.combined[['plates42_44']]['GMOY001603.T205', 'Ct1'] + plates.combined[['plates42_44']]['GMOY001603.T205', 'Ct3']) / 2
plates.combined[['plates42_44']]['GMOY001603.T205', 'Normalised.Ct'] <- (plates.combined[['plates42_44']]['GMOY001603.T205', 'Norm.Ct1'] + plates.combined[['plates42_44']]['GMOY001603.T205', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates53_55']]['GMOY001603.T239', 'Norm.Ct1'], plates.combined[['plates53_55']]['GMOY001603.T239', 'Norm.Ct2'], plate.147.group2['GMOY001603.T239', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates53_55']]$Ct3 <- plate.147.group2[rownames(plates.combined[['plates53_55']]), 'Ct']
plates.combined[['plates53_55']]$Norm.Ct3 <- plate.147.group2[rownames(plates.combined[['plates53_55']]), 'Normalised.Ct']
plates.combined[['plates53_55']]['GMOY001603.T239', 'Ct'] <- (plates.combined[['plates53_55']]['GMOY001603.T239', 'Ct1'] + plates.combined[['plates53_55']]['GMOY001603.T239', 'Ct3']) / 2
plates.combined[['plates53_55']]['GMOY001603.T239', 'Normalised.Ct'] <- (plates.combined[['plates53_55']]['GMOY001603.T239', 'Norm.Ct1'] + plates.combined[['plates53_55']]['GMOY001603.T239', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates90_92']]['GMOY001603.T369', 'Norm.Ct1'], plates.combined[['plates90_92']]['GMOY001603.T369', 'Norm.Ct2'], plate.147.group2['GMOY001603.T369', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates90_92']]$Ct3 <- plate.147.group2[rownames(plates.combined[['plates90_92']]), 'Ct']
plates.combined[['plates90_92']]$Norm.Ct3 <- plate.147.group2[rownames(plates.combined[['plates90_92']]), 'Normalised.Ct']
plates.combined[['plates90_92']]['GMOY001603.T369', 'Ct'] <- (plates.combined[['plates90_92']]['GMOY001603.T369', 'Ct1'] + plates.combined[['plates90_92']]['GMOY001603.T369', 'Ct3']) / 2
plates.combined[['plates90_92']]['GMOY001603.T369', 'Normalised.Ct'] <- (plates.combined[['plates90_92']]['GMOY001603.T369', 'Norm.Ct1'] + plates.combined[['plates90_92']]['GMOY001603.T369', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates105_107']]['GMOY001603.T418', 'Norm.Ct1'], plates.combined[['plates105_107']]['GMOY001603.T418', 'Norm.Ct2'], plate.147.group2['GMOY001603.T418', 'Normalised.Ct'])
# The third replicate is closest to the second, so we take the average of those two. 
plates.combined[['plates105_107']]$Ct3 <- plate.147.group2[rownames(plates.combined[['plates105_107']]), 'Ct']
plates.combined[['plates105_107']]$Norm.Ct3 <- plate.147.group2[rownames(plates.combined[['plates105_107']]), 'Normalised.Ct']
plates.combined[['plates105_107']]['GMOY001603.T418', 'Ct'] <- (plates.combined[['plates105_107']]['GMOY001603.T418', 'Ct2'] + plates.combined[['plates105_107']]['GMOY001603.T418', 'Ct3']) / 2
plates.combined[['plates105_107']]['GMOY001603.T418', 'Normalised.Ct'] <- (plates.combined[['plates105_107']]['GMOY001603.T418', 'Norm.Ct2'] + plates.combined[['plates105_107']]['GMOY001603.T418', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates137_139']]['GMOY001603.T534', 'Norm.Ct1'], plates.combined[['plates137_139']]['GMOY001603.T534', 'Norm.Ct2'], plate.147.group2['GMOY001603.T534', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates137_139']]$Ct3 <- plate.147.group2[rownames(plates.combined[['plates137_139']]), 'Ct']
plates.combined[['plates137_139']]$Norm.Ct3 <- plate.147.group2[rownames(plates.combined[['plates137_139']]), 'Normalised.Ct']
plates.combined[['plates137_139']]['GMOY001603.T534', 'Ct'] <- (plates.combined[['plates137_139']]['GMOY001603.T534', 'Ct1'] + plates.combined[['plates137_139']]['GMOY001603.T534', 'Ct3']) / 2
plates.combined[['plates137_139']]['GMOY001603.T534', 'Normalised.Ct'] <- (plates.combined[['plates137_139']]['GMOY001603.T534', 'Norm.Ct1'] + plates.combined[['plates137_139']]['GMOY001603.T534', 'Norm.Ct3']) / 2

# Now let's go to plate 148
plate.148 <- read.table('qPCR/191004/eric-tsetse-p148-191004 - Tabular Results.txt', sep = '\t', header = T, row.names = 1, na.strings = 'No Cq', quote = '', strip.white = T, col.names = c('Well', 'Ct', 'Tm', 'Gene'))
# Let's start with the repeats of plates 1-4 first
plate.148.genes <- c(rep(c('GMOY002920', hk.genes), 3),
                     c('GMOY003588', hk.genes),
                     rep(c('GMOY009908', hk.genes), 3),
                     c('GMOY009908', 'GMOY002920', 'GMOY003588')
                    )
plate.148$Gene <- plate.148.genes
s <- c(rep(c('T162', 'T236', 'T382', 'T456', 'T432', 'T502'), each = 3), 'T162', rep('NTC', 5))
plate.148$Sample <- factor(s)
plate.148 <- normalise.plate(plate.148)
rownames(plate.148) <- interaction(plate.148$Gene, plate.148$Sample)
# For each repeat, let's compare with the two originals
c(plates.combined[['plates30_32']]['GMOY002920.T162', 'Norm.Ct1'], plates.combined[['plates30_32']]['GMOY002920.T162', 'Norm.Ct2'], plate.148['GMOY002920.T162', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates30_32']]$Ct3 <- plate.148[rownames(plates.combined[['plates30_32']]), 'Ct']
plates.combined[['plates30_32']]$Norm.Ct3 <- plate.148[rownames(plates.combined[['plates30_32']]), 'Normalised.Ct']
plates.combined[['plates30_32']]['GMOY002920.T162', 'Ct'] <- (plates.combined[['plates30_32']]['GMOY002920.T162', 'Ct1'] + plates.combined[['plates30_32']]['GMOY002920.T162', 'Ct3']) / 2
plates.combined[['plates30_32']]['GMOY002920.T162', 'Normalised.Ct'] <- (plates.combined[['plates30_32']]['GMOY002920.T162', 'Norm.Ct1'] + plates.combined[['plates30_32']]['GMOY002920.T162', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates53_55']]['GMOY002920.T236', 'Norm.Ct1'], plates.combined[['plates53_55']]['GMOY002920.T236', 'Norm.Ct2'], plate.148['GMOY002920.T236', 'Normalised.Ct'])
# The third replicate is closest to the second, so we take the average of those two. 
plates.combined[['plates53_55']]$Ct3 <- plate.148[rownames(plates.combined[['plates53_55']]), 'Ct']
plates.combined[['plates53_55']]$Norm.Ct3 <- plate.148[rownames(plates.combined[['plates53_55']]), 'Normalised.Ct']
plates.combined[['plates53_55']]['GMOY002920.T236', 'Ct'] <- (plates.combined[['plates53_55']]['GMOY002920.T236', 'Ct2'] + plates.combined[['plates53_55']]['GMOY002920.T236', 'Ct3']) / 2
plates.combined[['plates53_55']]['GMOY002920.T236', 'Normalised.Ct'] <- (plates.combined[['plates53_55']]['GMOY002920.T236', 'Norm.Ct2'] + plates.combined[['plates53_55']]['GMOY002920.T236', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates94_96']]['GMOY002920.T382', 'Norm.Ct1'], plates.combined[['plates94_96']]['GMOY002920.T382', 'Norm.Ct2'], plate.148['GMOY002920.T382', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates94_96']]$Ct3 <- plate.148[rownames(plates.combined[['plates94_96']]), 'Ct']
plates.combined[['plates94_96']]$Norm.Ct3 <- plate.148[rownames(plates.combined[['plates94_96']]), 'Normalised.Ct']
plates.combined[['plates94_96']]['GMOY002920.T382', 'Ct'] <- (plates.combined[['plates94_96']]['GMOY002920.T382', 'Ct1'] + plates.combined[['plates94_96']]['GMOY002920.T382', 'Ct3']) / 2
plates.combined[['plates94_96']]['GMOY002920.T382', 'Normalised.Ct'] <- (plates.combined[['plates94_96']]['GMOY002920.T382', 'Norm.Ct1'] + plates.combined[['plates94_96']]['GMOY002920.T382', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates114_116']]['GMOY003588.T456', 'Norm.Ct1'], plates.combined[['plates114_116']]['GMOY003588.T456', 'Norm.Ct2'], plate.148['GMOY003588.T456', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates114_116']]$Ct3 <- plate.148[rownames(plates.combined[['plates114_116']]), 'Ct']
plates.combined[['plates114_116']]$Norm.Ct3 <- plate.148[rownames(plates.combined[['plates114_116']]), 'Normalised.Ct']
plates.combined[['plates114_116']]['GMOY003588.T456', 'Ct'] <- (plates.combined[['plates114_116']]['GMOY003588.T456', 'Ct1'] + plates.combined[['plates114_116']]['GMOY003588.T456', 'Ct3']) / 2
plates.combined[['plates114_116']]['GMOY003588.T456', 'Normalised.Ct'] <- (plates.combined[['plates114_116']]['GMOY003588.T456', 'Norm.Ct1'] + plates.combined[['plates114_116']]['GMOY003588.T456', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates109_111']]['GMOY009908.T432', 'Norm.Ct1'], plates.combined[['plates109_111']]['GMOY009908.T432', 'Norm.Ct2'], plate.148['GMOY009908.T432', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates109_111']]$Ct3 <- plate.148[rownames(plates.combined[['plates109_111']]), 'Ct']
plates.combined[['plates109_111']]$Norm.Ct3 <- plate.148[rownames(plates.combined[['plates109_111']]), 'Normalised.Ct']
plates.combined[['plates109_111']]['GMOY009908.T432', 'Ct'] <- (plates.combined[['plates109_111']]['GMOY009908.T432', 'Ct1'] + plates.combined[['plates109_111']]['GMOY009908.T432', 'Ct3']) / 2
plates.combined[['plates109_111']]['GMOY009908.T432', 'Normalised.Ct'] <- (plates.combined[['plates109_111']]['GMOY009908.T432', 'Norm.Ct1'] + plates.combined[['plates109_111']]['GMOY009908.T432', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates129_131']]['GMOY009908.T502', 'Norm.Ct1'], plates.combined[['plates129_131']]['GMOY009908.T502', 'Norm.Ct2'], plate.148['GMOY009908.T502', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates129_131']]$Ct3 <- plate.148[rownames(plates.combined[['plates129_131']]), 'Ct']
plates.combined[['plates129_131']]$Norm.Ct3 <- plate.148[rownames(plates.combined[['plates129_131']]), 'Normalised.Ct']
plates.combined[['plates129_131']]['GMOY009908.T502', 'Ct'] <- (plates.combined[['plates129_131']]['GMOY009908.T502', 'Ct1'] + plates.combined[['plates129_131']]['GMOY009908.T502', 'Ct3']) / 2
plates.combined[['plates129_131']]['GMOY009908.T502', 'Normalised.Ct'] <- (plates.combined[['plates129_131']]['GMOY009908.T502', 'Norm.Ct1'] + plates.combined[['plates129_131']]['GMOY009908.T502', 'Norm.Ct3']) / 2
#
c(plates.combined[['plates30_32']]['GMOY009908.T162', 'Norm.Ct1'], plates.combined[['plates30_32']]['GMOY009908.T162', 'Norm.Ct2'], plate.148['GMOY009908.T162', 'Normalised.Ct'])
# The third replicate is closest to the first, so we take the average of those two. 
plates.combined[['plates30_32']]$Ct3 <- plate.148[rownames(plates.combined[['plates30_32']]), 'Ct']
plates.combined[['plates30_32']]$Norm.Ct3 <- plate.148[rownames(plates.combined[['plates30_32']]), 'Normalised.Ct']
plates.combined[['plates30_32']]['GMOY009908.T162', 'Ct'] <- (plates.combined[['plates30_32']]['GMOY009908.T162', 'Ct1'] + plates.combined[['plates30_32']]['GMOY009908.T162', 'Ct3']) / 2
plates.combined[['plates30_32']]['GMOY009908.T162', 'Normalised.Ct'] <- (plates.combined[['plates30_32']]['GMOY009908.T162', 'Norm.Ct1'] + plates.combined[['plates30_32']]['GMOY009908.T162', 'Norm.Ct3']) / 2

# Now we have good data, let's output them to a csv file.

# Write a function that will turn a table where each row is a well to one where each row is a sample
change.table <- function(x){
	output.table <- data.frame(row.names = levels(x$Sample))
	for (g in genes)
		output.table[,g] <- tapply(subset(x, Gene == g)$Ct, subset(x, Gene == g)$Sample, function(x) x)[rownames(output.table)]
	for (g in age.genes)
		output.table[,paste(g, 'n', sep = '')] <- tapply(subset(x, Gene == g)$Normalised.Ct, subset(x, Gene == g)$Sample, function(x) x)[rownames(output.table)]
	output.table
}

for (p in names(plates.combined))
	plates[[p]] <- change.table(plates.combined[[p]])

all.plates <- do.call(rbind, plates)
all.plates <- all.plates[!grepl('NTC', rownames(all.plates)), ]
rownames(all.plates) <- sub('.*\\.', '', rownames(all.plates))

write.table(all.plates, file = 'tables/tsetse_qpcr.csv', sep = '\t', col.names = NA, quote = F)


