options(width=150)

root.path = function(...,create=FALSE){ # find project root above /code/ & build path from ...
  root = strsplit(getwd(),file.path('','code',''))[[1]][1]
  path = file.path(root,...)
  if (create){ suppressWarnings({ dir.create(path,recursive=TRUE) }) }
  return(path)
}
figname = function(name,...){ # e.g. /out/fig/fsa/.../name.pdf
  path = root.path('out','fig','fsa',...,create=TRUE)
  return(file.path(path,paste0(name,'.pdf')))
}
set.config = function(mode='10x10'){ # set config stuff in global list variable
  config = list( # see [[mode]] below
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
  config$c.type = c(
    'Household Contacts'     = 'home',
    'Non-Household Contacts' = 'other'
  )
    # add config stuff that doesn't depend on mode
  config = c(config,list(
    RR.C.global = c( # global scaling factor for contacts
      'home'  = 1,
      'other' = 1
    ),
    RR.C.decile = c(1 # TODO
      # 13.029573,9.067338,6.610766,5.409403,4.522559,
      #  4.068002,3.286063,2.687275,1.839832,1.000000
    ),
    phi = c( # odds of mobility vs observed devices
      'unobs.device' = .9,
      'no.device'    = .9
    ),
    h.y = c( # proportion of contacts formed at home
      'home'  = 1,
      'other' = 0
    ),
    t.ref   = c('2020-01','2020-02'), # REF period
    t.covid = c('2020-03','2020-04','2020-05','2020-06', # covid periods
                '2020-07','2020-08','2020-09','2020-10','2020-11','2020-12','2021-01'),
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
    ),
    labels = list( # convenience for plotting labels
      g = list(y='Home Decile (g)',x='Other Decile (g\')'),
      n = list(y='Home FSA (n)',   x='Other FSA (n\')'),
      a = list(y='Self Age (a)',   x='Other Age (a\')')
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
  config$age.max = max(config$age)+5
  config <<- config # make global
}
set.config()