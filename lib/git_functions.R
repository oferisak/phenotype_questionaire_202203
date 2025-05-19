library(usethis)

setup_git<-function(){
  ## now, working inside "the package", initialize git repository
  use_git()
}

create_git_repo<-function(){
  ## create github repository and configure as git remote
  use_github()
}

update_git_token<-function(){
  gitcreds::gitcreds_set()
}