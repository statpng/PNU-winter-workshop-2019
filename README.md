# kogo-2018
## KOGO HTS workshop website for 2018

See: https://htsworkshop.github.io/kogo-2018/

**Overview:** This is a simple [jekell based static site](http://jekyllrb.com/docs/home/). To view locally on your own machine (i.e. before pushing or submitting a pull 
request to this [htsworkshop GitHub](https://htsworkshop.github.io/kogo-2018/) repo) 
you will need to have the **jekyll** and **github-pages** gem setup (see further 
below for full instructions)


Consider updating RubyGems first (likely need sudo for these).

	sudo gem update --system

Then install the Jekyll Gem and the GitHub Gem

	gem install jekyll
	gem install github-pages

Optional: Pygments python based syntax highlighter

	pip install Pygments


## Basics of Jekyll websites
Jekyll websites are configured based on the contents of the various underscore prefixed files and folders. You can find out more about these here: http://jekyllrb.com/docs/structure/

However, most likely you will want to leave most of these alone and just add  
content to the day{2,3,4,5}.md files and create new files in the **class-material/** 
directory (i.e. add lecture slides, handouts, cheat-sheets etc.)

Please remember that all content is on the **master** branch! 
So you will want to be working on this branch and push back to this branch.

A typical workflow for folks that have been added as **"Collaborators"** would look something like this:

	## One time only clone
	git clone https://github.com/htsworkshop/kogo-2018.git
	cd web-2015

	## Edit your files (e.g. day2.md)
	vi day2.md

	## Check changes localy
	jekyll serve

	## Pull recent changes
	git pull origin master

	## Stage, commit and push your changes
	git status
	git add day2.md
	git commit -m "Your msg about changes"
	git push origin master


## How this site was built
Basic setup entailed:

	jekyll new 2015
	cd 2015
	## Edited site title, description etc. in _config.yml
	vi _config.yml  
	rm -rf _posts/   ##  we are not going to have a blog
	mv index.html blog.html  ## can delete later

Create a simple index.md file and have a quick look with 'jekyll serve'

	jekyll serve

After some more content addition I then followed the [instructions for adding 
to GitHub pages](http://jekyllrb.com/docs/github-pages/).

