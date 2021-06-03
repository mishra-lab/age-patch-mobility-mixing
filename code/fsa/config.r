suppressMessages({library(socialmixr)})
options(width=150)

DEBUG = TRUE
MODE = '10x10' # 10x10 = 10 deciles & 10+ age groups; 2x2 = 2 deciles & age groups
fig.ext = '.pdf'

root.path = function(...){
  root = strsplit(getwd(),file.path('','code',''))[[1]][1]
  return(file.path(root,...))
}
figname = function(name,...){
  path = root.path('out','fig','fsa',ifelse(DEBUG,'.debug',''),...)
  suppressWarnings({dir.create(path,recursive=TRUE)})
  return(file.path(path,paste0(name,fig.ext)))
}

info = list(
  '10x10' = list(
    age = c(
      '<12'   = 0,
      '12-16' = 12,
      '17-39' = 17,
      '40-44' = 40,
      '45-49' = 45,
      '50-54' = 50,
      '55-59' = 55,
      '60-64' = 60,
      '65-69' = 65,
      '70-74' = 70,
      '75-79' = 75,
      '80+'   = 80),
    c.type = c(
      'Household Contacts'     = 'Home',
      'Non-Household Contacts' = 'Other'),
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
      '16-59' = 16,
      '60+'   = 60),
    c.type = c(
      'Household Contacts'     = 'Home',
      'Non-Household Contacts' = 'Other'),
    group = c(
      '1-2'  = 1,
      '3-10' = 3)
))[[MODE]]
N = list(
  a = length(info$age),
  y = length(info$c.type),
  g = length(info$group)
)
eps = c(
  'home'  = 0.2354799,
  'other' = 0.1818460)
phi = c(
  'unobs.device' = .9,
  'no.device'    = .9
)
h.y = c(
  'home'  = 1,
  'other' = 0
)
mo.ref = c('2020-01','2020-02')
mo.covid = c('2020-03','2020-04','2020-05','2020-06',
  '2020-07','2020-08','2020-09','2020-10','2020-11','2020-12',
  '2021-01')
X.names = list(
  'g'  = names(info$group),
  'a'  = names(info$age),
  'g.' = names(info$group),
  'a.' = names(info$age)
)
labels = list(
  g = list(y='Home Decile (g)',x='Other Decile (g\')'),
  n = list(y='Home FSA (n)',   x='Other FSA (n\')'),
  a = list(y='Self Age (a)',   x='Other Age (a\')')
)
age.contact = c(
  '00-04' =  0 , '05-09' =  5,
  '10-14' = 10 , '15-19' = 15,
  '20-24' = 20 , '25-29' = 25,
  '30-34' = 30 , '35-39' = 35,
  '40-44' = 40 , '45-49' = 45,
  '50-54' = 50 , '55-59' = 55,
  '60-64' = 60 , '65-69' = 65,
  '70-74' = 70 , '75-80' = 75) # TODO: double check defs
