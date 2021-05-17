suppressMessages({library(socialmixr)})
options(width=150)

fig.ext = '.pdf'

root.path = function(...){
  root = strsplit(getwd(),file.path('','code',''))[[1]][1]
  return(file.path(root,...))
}
figname = function(name,...){
  root.path('out','fig','fsa',...,paste0(name,fig.ext))
}

info = list(
  age = c(
    '<12'   = 0,
    '13-16' = 13,
    '17-39' = 17,
    '40-44' = 40,
    '45-49' = 45,
    '50-54' = 50,
    '55-59' = 55,
    '60-64' = 60,
    '65-69' = 65,
    '70-74' = 70,
    '75-79' = 75,
    '80+  ' = 80),
  c.type = c(
    'Home'  = 'Home',
    'Other' = 'Other'),
  decile = c(1:10)
)
names(info$decile) = info$decile
N = list(
  a = length(info$age),
  y = length(info$c.type),
  g = length(info$decile)
)
eps = c(
  'home'  = 0.08785141,
  'other' = 0.18124466
)
OR.travel.unobs = .5
h.y = c(
  'home'  = 1,
  'other' = 0
)
mo.ref = c('2020-01','2020-02')
X.names = list(
  'g'  = names(info$decile),
  'a'  = names(info$age),
  'g.' = names(info$decile),
  'a.' = names(info$age)
)
labels = list(
  g = list(y='Index Decile (g)',x='Other Decile (g\')'),
  a = list(y='Index Age (a)',   x='Other Age (a\')'),
  n = list(y='Index FSA (n)',   x='Other FSA (n\')')
)
