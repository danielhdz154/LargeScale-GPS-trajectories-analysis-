--DROP TABLE template_sygic;

-- Create schema

CREATE TABLE template_sygic (
  log_time timestamp,
  latitude real,
  longitude real,
  x1 int,
  speed real,
  X2 int,
  x3 int ,
  x4 int,
  x5 int,
  speed_2 real,
  x6 int,
  x7 int,
  region_code varchar(5),
  session_id varchar(100),
  device_id varchar(100),
  platform varchar(10),
  device_m varchar(100),
  app_v varchar(100),
  unknown_id varchar(100)
);

-- import from CSV  # CHANGE PATH AND CSV FILE NAME
copy template_sygic from 'C:\Users\Daniel Hernandez\Documents\UNIVERSITY OF MELBOURNE\Subjects\Research Project\Data\Sygic\SygicAustralia2016-04-29.csv' DELIMITER ',' CSV encoding 'UTF8';

-- Remove unecessary columns and change name of table
ALTER TABLE template_sygic
    drop column x1, drop column x2, DROP column x3, DROP column x4, DROP column x5, DROP column x6, drop column x7,
    drop column app_v, drop column device_m, drop column platform;

ALTER TABLE template_sygic RENAME TO sygic_melb_2016_04_29;

-- Extract points for Melbourne

DELETE FROM sygic_melb_2016_04_24 WHERE latitude<(-38.5072) or latitude>(-37.4209) or longitude<144.5660 or longitude>145.5054;

-- Add geometry and convert to UTM coordinates

SELECT AddGeometryColumn('public', 'sygic_melb_2016_04_24', 'geom', 28355, 'POINT', 2, false)
UPDATE sygic_melb_2016_04_24 SET geom = ST_Transform(ST_SetSRID(ST_MakePoint(longitude, latitude), 4326),28355);

-- Create INDEX for session_id #CHANGE NAME traj_index_04_29
CREATE INDEX traj_index_04_29 ON sygic_melb_2016_04_29(session_id DESC NULLS LAST);
