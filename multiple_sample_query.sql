select tissue_info.*, demographics.species, demographics.breed, path.diagnosis, path.infectious from tissue_info
	INNER JOIN (SELECT ahl_id
    from tissue_info
    GROUP BY ahl_id
    HAVING COUNT(type)>3) dup
    on tissue_info.ahl_id = dup.ahl_id
    JOIN demographics on tissue_info.ahl_id=demographics.ahl_id
    JOIN path on tissue_info.ahl_id=path.ahl_id
		where (demographics.species='Horse' or demographics.species='horse') and (type='RNA')
		into outfile '~/test.csv'
		fields terminated by ',';
