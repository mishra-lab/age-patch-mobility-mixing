suppressMessages({library(socialmixr)})
options(width=150)

root.path = function(...){
  root = strsplit(getwd(),file.path('','code',''))[[1]][1]
  return(file.path(root,...))
}

info = list(
  age = c(
    '<16'   = 0,
    '16-39' = 16,
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
    'Home'      = 'home',
    'Work'      = 'work',
    'School'    = 'school',
    'Transport' = 'transport',
    'Leisure'   = 'leisure',
    'Other'     = 'otherplace'),
  decile = c(1:10)
)
names(info$decile) = info$decile
N = list(
  a = length(info$age),
  y = length(info$c.type),
  g = length(info$decile)
)
eps = c(
  'home'       = 0.08805234,
  'work'       = 0.02565621,
  'school'     = 0.22177813,
  'transport'  = 0.12430744,
  'leisure'    = 0.13311496,
  'otherplace' = 0.05587195
)
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
