# install.packages('rsconnect')

# https://www.shinyapps.io/admin/#/dashboard
rsconnect::setAccountInfo(name=' ',
                          token=' ',
                          secret='<SECRET>')

# deploy your app
rsconnect::deployApp(
  appDir = ".",
  appPrimaryDoc = "app.R"
)
