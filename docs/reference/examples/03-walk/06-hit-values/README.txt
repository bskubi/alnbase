How a hit row's values are taken from columns.
flat.toml: anchor on column 0 (default) vs column 1 (mark ".+").
lists.toml: extra captures ('^') switch the WHOLE file to list columns;
  captures with no read base get NULL qual/off_5p/off_3p.
tags.toml: BAM tag then `alnbase extract`; extract fills refr_base only
  when the query pins the anchor's reference symbol (NULL for '~').
