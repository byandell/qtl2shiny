# README #

Yandell R/qtl2shiny project.

### What is this repository for? ###

* Code to create shiny interface for [R/doqtl2](https://bitbucket.org/byandell/doqtl2).
* Version 0.2

### What has been done ###

* Shiny interface for multiple phenotypes with modules and dashboard
  + invoke shiny dashboard as doqtl2_app()
  + download RDS data
  + create plots and scans on the fly
* organized in package R/qtl2shiny
* depends on R/doqtl2

### What are open development issues ###

Major issues

* Navigation across panels imporved
* keeping dependencies in reactives straight
* speed up access
* save intermediate calculations that are reused

Minor issues

* shiny
  + reveal plots more (settings on sidebar? tabs?)
  + user save settings for quick replay of shiny
* markdown
  + flexdashboard for Rmd to create dynamic reports
  + link Rmd and shiny

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* [Brian Yandell](http://bitbucket.org/byandell)
* [Karl Broman](http://bitbucket.org/kbroman)
  + [R/qtl2](http://kbroman.org/qtl2/)
