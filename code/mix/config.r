options(width=150)

root.path = function(...,create=FALSE){ # find project root above /code/ & build path from ...
  root = strsplit(getwd(),file.path('','code',''))[[1]][1]
  path = file.path(root,...)
  if (create){ suppressWarnings({ dir.create(path,recursive=TRUE) }) }
  return(path)
}
figname = function(name,...,ext='.pdf'){ # e.g. /out/fig/fsa/.../name.pdf
  path = root.path('out','fig','mix',...,create=TRUE)
  return(file.path(path,paste0(name,ext)))
}
set.config = function(mode='10x10',n.y=2,...){ # set config stuff in global list variable
  config = list( # see [[mode]] below
    # TODO: remove support for 2x2 ?
    '10x10' = list(
      age = c(
        '<12'   = 0,
        '12-15' = 12,
        '16-39' = 16,
        '40-44' = 40,
        '45-49' = 45,
        '50-54' = 50,
        '55-59' = 55,
        '60-64' = 60,
        '65-69' = 65,
        '70-74' = 70,
        '75-79' = 75,
        '80+'   = 80),
      group = c(
        '1'  = 1,
        '2'  = 2,
        '3'  = 3,
        '4'  = 4,
        '5'  = 5,
        '6'  = 6,
        '7'  = 7,
        '8'  = 8,
        '9'  = 9,
        '10' = 10)
    ),
    '2x2' = list(
      age = c(
        '12-15' = 12,
        '16-59' = 16,
        '60+'   = 60),
      group = c(
        '1-2'  = 1,
        '3-10' = 3)
  ))[[mode]]
  config$mode = mode
  config = c(config,list(
    '2' = list(
      c.map = list(hh='home',nhh=c('work','school','others')),
      c.type = c(
        'Household Contacts'     = 'hh',
        'Non-Household Contacts' = 'nhh'
      ),
      h.y       = c('hh'=1,'nhh'=0),
      RC.global = c('hh'=1,'nhh'=0)
    ),
    '4' = list(
      c.map=list(home='home',work='work',school='school',others='others'),
      c.type = c(
        'Home'   = 'home',
        'Work'   = 'work',
        'School' = 'school',
        'Other'  = 'others'
      ),
      h.y       = c('home'=1,'work'=0,'school'=0,'others'=0),
      RC.global = c('home'=1,'work'=1,'school'=1,'others'=1)
    )
  )[[paste(n.y)]])
  # add config stuff that doesn't depend on mode
  config = c(config,list(
    method = list(
      'home.pool'   = TRUE,
      'age.adapt'   = TRUE,
      'age.by.type' = TRUE
    ),
    RC.decile = c(1
      # 7.643163, 5.755709, 4.448017, 3.846107, 3.246633,
      # 2.983167, 2.539929, 2.136271, 1.625427, 1.000000
      # 13.029573, 9.067338, 6.610766, 5.409403, 4.522559,
      #  4.068002, 3.286063, 2.687275, 1.839832, 1.000000
    ),
    phi = c( # odds of mobility vs observed devices
      'unobs.device' = .9,
      'no.device'    = .9
    ),
    t.ref   = c('2020-01','2020-02'), # REF period
    t.covid = c('2020-03','2020-04','2020-05','2020-06', # covid periods
                '2020-07','2020-08','2020-09','2020-10','2020-11','2020-12'),
    labels = list( # convenience for plotting labels
      g = list(x='Home Decile (g)',y='Other Decile (g\')'),
      n = list(x='Home FSA (n)',   y='Other FSA (n\')'),
      a = list(x='Self Age (a)',   y='Other Age (a\')')
    ),
    age.contact = c( # Prem 2017 POLYMOD age groups
      '00-04' =  0 , '05-09' =  5,
      '10-14' = 10 , '15-19' = 15,
      '20-24' = 20 , '25-29' = 25,
      '30-34' = 30 , '35-39' = 35,
      '40-44' = 40 , '45-49' = 45,
      '50-54' = 50 , '55-59' = 55,
      '60-64' = 60 , '65-69' = 65,
      '70-74' = 70 , '75-80' = 75)
  ))
  args = list(...)
  for (name in names(args)){
    config[[name]] = args[[name]]
  }
  config = c(config,list(
    X.names = list( # convenience list for naming matrix dims
      'g'  = names(config$group),
      'a'  = names(config$age),
      'g.' = names(config$group),
      'a.' = names(config$age)
    ),
    N = list( # convenience
      a = length(config$age),
      y = length(config$c.type),
      g = length(config$group)
    )
  ))
  config <<- config # make global
}
set.config()