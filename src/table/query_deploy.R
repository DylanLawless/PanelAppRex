# install.packages('rsconnect')

# https://www.shinyapps.io/admin/#/dashboard
rsconnect::setAccountInfo(name='switzerlandomics',
                          token='F852B8DEAD2012EBD9BB61D3AA4412A3',
                          secret='<SECRET>')

# deploy your app
rsconnect::deployApp(
  appDir = ".",
  appPrimaryDoc = "app.R"
)
