




<!DOCTYPE html>
<html class="   ">
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# object: http://ogp.me/ns/object# article: http://ogp.me/ns/article# profile: http://ogp.me/ns/profile#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    
    
    <title>MosaicHunter/README.md at master · AugustHuang/MosaicHunter</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <meta property="fb:app_id" content="1401488693436528"/>

      <meta content="@github" name="twitter:site" /><meta content="summary" name="twitter:card" /><meta content="AugustHuang/MosaicHunter" name="twitter:title" /><meta content="MosaicHunter - A tool for detecting postzygotic single-nucleotide mutations in human whole-genome sequencing data." name="twitter:description" /><meta content="https://avatars2.githubusercontent.com/u/8017324?s=400" name="twitter:image:src" />
<meta content="GitHub" property="og:site_name" /><meta content="object" property="og:type" /><meta content="https://avatars2.githubusercontent.com/u/8017324?s=400" property="og:image" /><meta content="AugustHuang/MosaicHunter" property="og:title" /><meta content="https://github.com/AugustHuang/MosaicHunter" property="og:url" /><meta content="MosaicHunter - A tool for detecting postzygotic single-nucleotide mutations in human whole-genome sequencing data." property="og:description" />

    <link rel="assets" href="https://assets-cdn.github.com/">
    <link rel="conduit-xhr" href="https://ghconduit.com:25035">
    <link rel="xhr-socket" href="/_sockets" />

    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
      <meta name="google-analytics" content="UA-3769691-2">

    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="collector-cdn.github.com" name="octolytics-script-host" /><meta content="github" name="octolytics-app-id" /><meta content="CA267E12:1805:494032:53B51566" name="octolytics-dimension-request_id" /><meta content="8017324" name="octolytics-actor-id" /><meta content="AugustHuang" name="octolytics-actor-login" /><meta content="f6e9ce0826cf24360f74f473a591c5869353daceae1ef97fd920d27f45e4d2fa" name="octolytics-actor-hash" />
    

    
    
    <link rel="icon" type="image/x-icon" href="https://assets-cdn.github.com/favicon.ico" />


    <meta content="authenticity_token" name="csrf-param" />
<meta content="EPVfP0crvTNSjUh7T1XLR0aMX/oX6ht7kznrcuPn4a0HeJ1YmqBqDUx0ynZYm6DcYOzgzl8fYHehhwt86FcldQ==" name="csrf-token" />

    <link href="https://assets-cdn.github.com/assets/github-a3943029fb2330481c4a6367eccd68e84b5cb8d7.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://assets-cdn.github.com/assets/github2-36643cb205883fc24ac19c4277fa9a1e2b53182b.css" media="all" rel="stylesheet" type="text/css" />
    


    <meta http-equiv="x-pjax-version" content="3848826c0acf5ea8e33de67266b149ac">

      
  <meta name="description" content="MosaicHunter - A tool for detecting postzygotic single-nucleotide mutations in human whole-genome sequencing data." />


  <meta content="8017324" name="octolytics-dimension-user_id" /><meta content="AugustHuang" name="octolytics-dimension-user_login" /><meta content="21319407" name="octolytics-dimension-repository_id" /><meta content="AugustHuang/MosaicHunter" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="false" name="octolytics-dimension-repository_is_fork" /><meta content="21319407" name="octolytics-dimension-repository_network_root_id" /><meta content="AugustHuang/MosaicHunter" name="octolytics-dimension-repository_network_root_nwo" />

  <link href="https://github.com/AugustHuang/MosaicHunter/commits/master.atom" rel="alternate" title="Recent Commits to MosaicHunter:master" type="application/atom+xml" />

  </head>


  <body class="logged_in  env-production windows vis-public page-blob">
    <a href="#start-of-content" tabindex="1" class="accessibility-aid js-skip-to-content">Skip to content</a>
    <div class="wrapper">
      
      
      
      


      <div class="header header-logged-in true">
  <div class="container clearfix">

    <a class="header-logo-invertocat" href="https://github.com/" aria-label="Homepage">
  <span class="mega-octicon octicon-mark-github"></span>
</a>


    
    <a href="/notifications" aria-label="You have no unread notifications" class="notification-indicator tooltipped tooltipped-s" data-hotkey="g n">
        <span class="mail-status all-read"></span>
</a>

      <div class="command-bar js-command-bar  in-repository">
          <form accept-charset="UTF-8" action="/search" class="command-bar-form" id="top_search_form" method="get">

<div class="commandbar">
  <span class="message"></span>
  <input type="text" data-hotkey="s, /" name="q" id="js-command-bar-field" placeholder="Search or type a command" tabindex="1" autocapitalize="off"
    
    data-username="AugustHuang"
      data-repo="AugustHuang/MosaicHunter"
      data-branch="master"
      data-sha="62b9f992982ac7ee4b505b59b9f8c77123448b50"
  >
  <div class="display hidden"></div>
</div>

    <input type="hidden" name="nwo" value="AugustHuang/MosaicHunter" />

    <div class="select-menu js-menu-container js-select-menu search-context-select-menu">
      <span class="minibutton select-menu-button js-menu-target" role="button" aria-haspopup="true">
        <span class="js-select-button">This repository</span>
      </span>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container" aria-hidden="true">
        <div class="select-menu-modal">

          <div class="select-menu-item js-navigation-item js-this-repository-navigation-item selected">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" class="js-search-this-repository" name="search_target" value="repository" checked="checked" />
            <div class="select-menu-item-text js-select-button-text">This repository</div>
          </div> <!-- /.select-menu-item -->

          <div class="select-menu-item js-navigation-item js-all-repositories-navigation-item">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" name="search_target" value="global" />
            <div class="select-menu-item-text js-select-button-text">All repositories</div>
          </div> <!-- /.select-menu-item -->

        </div>
      </div>
    </div>

  <span class="help tooltipped tooltipped-s" aria-label="Show command bar help">
    <span class="octicon octicon-question"></span>
  </span>


  <input type="hidden" name="ref" value="cmdform">

</form>
        <ul class="top-nav">
          <li class="explore"><a href="/explore">Explore</a></li>
            <li><a href="https://gist.github.com">Gist</a></li>
            <li><a href="/blog">Blog</a></li>
          <li><a href="https://help.github.com">Help</a></li>
        </ul>
      </div>

    


  <ul id="user-links">
    <li>
      <a href="/AugustHuang" class="name">
        <img alt="AugustHuang" class=" js-avatar" data-user="8017324" height="20" src="https://avatars1.githubusercontent.com/u/8017324?s=140" width="20" /> AugustHuang
      </a>
    </li>

    <li class="new-menu dropdown-toggle js-menu-container">
      <a href="#" class="js-menu-target tooltipped tooltipped-s" aria-label="Create new...">
        <span class="octicon octicon-plus"></span>
        <span class="dropdown-arrow"></span>
      </a>

      <div class="new-menu-content js-menu-content">
      </div>
    </li>

    <li>
      <a href="/settings/profile" id="account_settings"
        class="tooltipped tooltipped-s"
        aria-label="Account settings ">
        <span class="octicon octicon-tools"></span>
      </a>
    </li>
    <li>
      <form class="logout-form" action="/logout" method="post">
        <button class="sign-out-button tooltipped tooltipped-s" aria-label="Sign out">
          <span class="octicon octicon-sign-out"></span>
        </button>
      </form>
    </li>

  </ul>

<div class="js-new-dropdown-contents hidden">
  

<ul class="dropdown-menu">
  <li>
    <a href="/new"><span class="octicon octicon-repo"></span> New repository</a>
  </li>
  <li>
    <a href="/organizations/new"><span class="octicon octicon-organization"></span> New organization</a>
  </li>


    <li class="section-title">
      <span title="AugustHuang/MosaicHunter">This repository</span>
    </li>
      <li>
        <a href="/AugustHuang/MosaicHunter/issues/new"><span class="octicon octicon-issue-opened"></span> New issue</a>
      </li>
      <li>
        <a href="/AugustHuang/MosaicHunter/settings/collaboration"><span class="octicon octicon-person"></span> New collaborator</a>
      </li>
</ul>

</div>


    
  </div>
</div>

      

        



      <div id="start-of-content" class="accessibility-aid"></div>
          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    <div id="js-flash-container">
      
    </div>
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        
<ul class="pagehead-actions">

    <li class="subscription">
      <form accept-charset="UTF-8" action="/notifications/subscribe" class="js-social-container" data-autosubmit="true" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="sv1PKXdZNwhqmVeXnO8hs3HDLGjffmsanh73luYzUOg1qVSeS4iKpqTLyspsPz6DnW0Zm3kPniwhX/+r3EHSeA==" /></div>  <input id="repository_id" name="repository_id" type="hidden" value="21319407" />

    <div class="select-menu js-menu-container js-select-menu">
      <a class="social-count js-social-count" href="/AugustHuang/MosaicHunter/watchers">
        3
      </a>
      <span class="minibutton select-menu-button with-count js-menu-target" role="button" tabindex="0" aria-haspopup="true">
        <span class="js-select-button">
          <span class="octicon octicon-eye"></span>
          Unwatch
        </span>
      </span>

      <div class="select-menu-modal-holder">
        <div class="select-menu-modal subscription-menu-modal js-menu-content" aria-hidden="true">
          <div class="select-menu-header">
            <span class="select-menu-title">Notification status</span>
            <span class="octicon octicon-x js-menu-close"></span>
          </div> <!-- /.select-menu-header -->

          <div class="select-menu-list js-navigation-container" role="menu">

            <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_included" name="do" type="radio" value="included" />
                <h4>Not watching</h4>
                <span class="description">You only receive notifications for conversations in which you participate or are @mentioned.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye"></span>
                  Watch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item selected" role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input checked="checked" id="do_subscribed" name="do" type="radio" value="subscribed" />
                <h4>Watching</h4>
                <span class="description">You receive notifications for all conversations in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye"></span>
                  Unwatch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_ignore" name="do" type="radio" value="ignore" />
                <h4>Ignoring</h4>
                <span class="description">You do not receive any notifications for conversations in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-mute"></span>
                  Stop ignoring
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

          </div> <!-- /.select-menu-list -->

        </div> <!-- /.select-menu-modal -->
      </div> <!-- /.select-menu-modal-holder -->
    </div> <!-- /.select-menu -->

</form>
    </li>

  <li>
    

  <div class="js-toggler-container js-social-container starring-container ">

    <form accept-charset="UTF-8" action="/AugustHuang/MosaicHunter/unstar" class="js-toggler-form starred" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="CtAEkjnWYfRGzQ621u4l8Ykw2dq/APtKP0hyEiiodYhZnnCGTX+gW8GaKPKsuvWBg+OyLbwjZChdmUsQIDXSeQ==" /></div>
      <button
        class="minibutton with-count js-toggler-target star-button"
        aria-label="Unstar this repository" title="Unstar AugustHuang/MosaicHunter">
        <span class="octicon octicon-star"></span>
        Unstar
      </button>
        <a class="social-count js-social-count" href="/AugustHuang/MosaicHunter/stargazers">
          0
        </a>
</form>
    <form accept-charset="UTF-8" action="/AugustHuang/MosaicHunter/star" class="js-toggler-form unstarred" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="JUcrYBHAFFvU/dKUcCnFLmlGSm2I4v9QubSL2cb5SIQ0+XjSvzNGjXJtLZBR3r6IFVJjnVxyZ4epzF6LBF8kvQ==" /></div>
      <button
        class="minibutton with-count js-toggler-target star-button"
        aria-label="Star this repository" title="Star AugustHuang/MosaicHunter">
        <span class="octicon octicon-star"></span>
        Star
      </button>
        <a class="social-count js-social-count" href="/AugustHuang/MosaicHunter/stargazers">
          0
        </a>
</form>  </div>

  </li>


        <li>
          <a href="/AugustHuang/MosaicHunter/fork" class="minibutton with-count js-toggler-target fork-button lighter tooltipped-n" title="Fork your own copy of AugustHuang/MosaicHunter to your account" aria-label="Fork your own copy of AugustHuang/MosaicHunter to your account" rel="nofollow" data-method="post">
            <span class="octicon octicon-repo-forked"></span>
            Fork
          </a>
          <a href="/AugustHuang/MosaicHunter/network" class="social-count">0</a>
        </li>

</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="repo-label"><span>public</span></span>
          <span class="mega-octicon octicon-repo"></span>
          <span class="author"><a href="/AugustHuang" class="url fn" itemprop="url" rel="author"><span itemprop="title">AugustHuang</span></a></span><!--
       --><span class="path-divider">/</span><!--
       --><strong><a href="/AugustHuang/MosaicHunter" class="js-current-repository js-repo-home-link">MosaicHunter</a></strong>

          <span class="page-context-loader">
            <img alt="" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
          </span>

        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">
      <div class="repository-with-sidebar repo-container new-discussion-timeline js-new-discussion-timeline  ">
        <div class="repository-sidebar clearfix">
            

<div class="sunken-menu vertical-right repo-nav js-repo-nav js-repository-container-pjax js-octicon-loaders">
  <div class="sunken-menu-contents">
    <ul class="sunken-menu-group">
      <li class="tooltipped tooltipped-w" aria-label="Code">
        <a href="/AugustHuang/MosaicHunter" aria-label="Code" class="selected js-selected-navigation-item sunken-menu-item" data-hotkey="g c" data-pjax="true" data-selected-links="repo_source repo_downloads repo_commits repo_releases repo_tags repo_branches /AugustHuang/MosaicHunter">
          <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

        <li class="tooltipped tooltipped-w" aria-label="Issues">
          <a href="/AugustHuang/MosaicHunter/issues" aria-label="Issues" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-hotkey="g i" data-selected-links="repo_issues /AugustHuang/MosaicHunter/issues">
            <span class="octicon octicon-issue-opened"></span> <span class="full-word">Issues</span>
            <span class='counter'>0</span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>

      <li class="tooltipped tooltipped-w" aria-label="Pull Requests">
        <a href="/AugustHuang/MosaicHunter/pulls" aria-label="Pull Requests" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-hotkey="g p" data-selected-links="repo_pulls /AugustHuang/MosaicHunter/pulls">
            <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull Requests</span>
            <span class='counter'>0</span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


        <li class="tooltipped tooltipped-w" aria-label="Wiki">
          <a href="/AugustHuang/MosaicHunter/wiki" aria-label="Wiki" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-hotkey="g w" data-selected-links="repo_wiki /AugustHuang/MosaicHunter/wiki">
            <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
    </ul>
    <div class="sunken-menu-separator"></div>
    <ul class="sunken-menu-group">

      <li class="tooltipped tooltipped-w" aria-label="Pulse">
        <a href="/AugustHuang/MosaicHunter/pulse" aria-label="Pulse" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="pulse /AugustHuang/MosaicHunter/pulse">
          <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped tooltipped-w" aria-label="Graphs">
        <a href="/AugustHuang/MosaicHunter/graphs" aria-label="Graphs" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_graphs repo_contributors /AugustHuang/MosaicHunter/graphs">
          <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped tooltipped-w" aria-label="Network">
        <a href="/AugustHuang/MosaicHunter/network" aria-label="Network" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-selected-links="repo_network /AugustHuang/MosaicHunter/network">
          <span class="octicon octicon-repo-forked"></span> <span class="full-word">Network</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>
    </ul>


      <div class="sunken-menu-separator"></div>
      <ul class="sunken-menu-group">
        <li class="tooltipped tooltipped-w" aria-label="Settings">
          <a href="/AugustHuang/MosaicHunter/settings" aria-label="Settings" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_settings /AugustHuang/MosaicHunter/settings">
            <span class="octicon octicon-tools"></span> <span class="full-word">Settings</span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
      </ul>
  </div>
</div>

              <div class="only-with-full-nav">
                

  

<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=push">
  <h3><strong>HTTPS</strong> clone URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/AugustHuang/MosaicHunter.git" readonly="readonly">
    <span class="url-box-clippy">
    <button aria-label="Copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="https://github.com/AugustHuang/MosaicHunter.git" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>

  

<div class="clone-url "
  data-protocol-type="ssh"
  data-url="/users/set_protocol?protocol_selector=ssh&amp;protocol_type=push">
  <h3><strong>SSH</strong> clone URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="git@github.com:AugustHuang/MosaicHunter.git" readonly="readonly">
    <span class="url-box-clippy">
    <button aria-label="Copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="git@github.com:AugustHuang/MosaicHunter.git" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>

  

<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=push">
  <h3><strong>Subversion</strong> checkout URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/AugustHuang/MosaicHunter" readonly="readonly">
    <span class="url-box-clippy">
    <button aria-label="Copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="https://github.com/AugustHuang/MosaicHunter" data-copied-hint="Copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>


<p class="clone-options">You can clone with
      <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>,
      <a href="#" class="js-clone-selector" data-protocol="ssh">SSH</a>,
      or <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>.
  <a href="https://help.github.com/articles/which-remote-url-should-i-use" class="help tooltipped tooltipped-n" aria-label="Get help on which URL is right for you.">
    <span class="octicon octicon-question"></span>
  </a>
</p>


  <a href="http://windows.github.com" class="minibutton sidebar-button" title="Save AugustHuang/MosaicHunter to your computer and use it in GitHub Desktop." aria-label="Save AugustHuang/MosaicHunter to your computer and use it in GitHub Desktop.">
    <span class="octicon octicon-device-desktop"></span>
    Clone in Desktop
  </a>

                <a href="/AugustHuang/MosaicHunter/archive/master.zip"
                   class="minibutton sidebar-button"
                   aria-label="Download AugustHuang/MosaicHunter as a zip file"
                   title="Download AugustHuang/MosaicHunter as a zip file"
                   rel="nofollow">
                  <span class="octicon octicon-cloud-download"></span>
                  Download ZIP
                </a>
              </div>
        </div><!-- /.repository-sidebar -->

        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
          


<a href="/AugustHuang/MosaicHunter/blob/a4a73a24c100ca39d058939cc796911dbdd4f9e5/README.md" class="hidden js-permalink-shortcut" data-hotkey="y">Permalink</a>

<!-- blob contrib key: blob_contributors:v21:6ea8636f855127708747ba98d93860b0 -->

<p title="This is a placeholder element" class="js-history-link-replace hidden"></p>

<div class="file-navigation">
  

<div class="select-menu js-menu-container js-select-menu" >
  <span class="minibutton select-menu-button js-menu-target css-truncate" data-hotkey="w"
    data-master-branch="master"
    data-ref="master"
    title="master"
    role="button" aria-label="Switch branches or tags" tabindex="0" aria-haspopup="true">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button css-truncate-target">master</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax aria-hidden="true">

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-x js-menu-close"></span>
      </div> <!-- /.select-menu-header -->

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Find or create a branch…" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Find or create a branch…">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div><!-- /.select-menu-tabs -->
      </div><!-- /.select-menu-filters -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/AugustHuang/MosaicHunter/blob/master/README.md"
                 data-name="master"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text css-truncate-target"
                 title="master">master</a>
            </div> <!-- /.select-menu-item -->
        </div>

          <form accept-charset="UTF-8" action="/AugustHuang/MosaicHunter/branches" class="js-create-branch select-menu-item select-menu-new-item-form js-navigation-item js-new-item-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="nWx3gZeZPuw48jOEhHDgsTeR2CByoOYB2ztulGK2FeDjZ1XQt6qaEQ30+ayDOF/ukdh0WOo7utBsodN4jNucig==" /></div>
            <span class="octicon octicon-git-branch select-menu-item-icon"></span>
            <div class="select-menu-item-text">
              <h4>Create branch: <span class="js-new-item-name"></span></h4>
              <span class="description">from ‘master’</span>
            </div>
            <input type="hidden" name="name" id="name" class="js-new-item-value">
            <input type="hidden" name="branch" id="branch" value="master" />
            <input type="hidden" name="path" id="path" value="README.md" />
          </form> <!-- /.select-menu-item -->

      </div> <!-- /.select-menu-list -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

    </div> <!-- /.select-menu-modal -->
  </div> <!-- /.select-menu-modal-holder -->
</div> <!-- /.select-menu -->

  <div class="button-group right">
    <a href="/AugustHuang/MosaicHunter/find/master"
          class="js-show-file-finder minibutton empty-icon tooltipped tooltipped-s"
          data-pjax
          data-hotkey="t"
          aria-label="Quickly jump between files">
      <span class="octicon octicon-list-unordered"></span>
    </a>
    <button class="js-zeroclipboard minibutton zeroclipboard-button"
          data-clipboard-text="README.md"
          aria-label="Copy to clipboard"
          data-copied-hint="Copied!">
      <span class="octicon octicon-clippy"></span>
    </button>
  </div>

  <div class="breadcrumb">
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/AugustHuang/MosaicHunter" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">MosaicHunter</span></a></span></span><span class="separator"> / </span><strong class="final-path">README.md</strong>
  </div>
</div>


  <div class="commit file-history-tease">
      <img alt="AugustHuang" class="main-avatar js-avatar" data-user="8017324" height="24" src="https://avatars1.githubusercontent.com/u/8017324?s=140" width="24" />
      <span class="author"><a href="/AugustHuang" rel="author">AugustHuang</a></span>
      <time datetime="2014-07-03T14:55:29+08:00" is="relative-time">July 03, 2014</time>
      <div class="commit-title">
          <a href="/AugustHuang/MosaicHunter/commit/a4a73a24c100ca39d058939cc796911dbdd4f9e5" class="message" data-pjax="true" title="Revised README">Revised README</a>
      </div>

    <div class="participation">
      <p class="quickstat"><a href="#blob_contributors_box" rel="facebox"><strong>1</strong>  contributor</a></p>
      
    </div>
    <div id="blob_contributors_box" style="display:none">
      <h2 class="facebox-header">Users who have contributed to this file</h2>
      <ul class="facebox-user-list">
          <li class="facebox-user-list-item">
            <img alt="AugustHuang" class=" js-avatar" data-user="8017324" height="24" src="https://avatars1.githubusercontent.com/u/8017324?s=140" width="24" />
            <a href="/AugustHuang">AugustHuang</a>
          </li>
      </ul>
    </div>
  </div>

<div class="file-box">
  <div class="file">
    <div class="meta clearfix">
      <div class="info file-name">
        <span class="icon"><b class="octicon octicon-file-text"></b></span>
        <span class="mode" title="File Mode">file</span>
        <span class="meta-divider"></span>
          <span>84 lines (71 sloc)</span>
          <span class="meta-divider"></span>
        <span>3.833 kb</span>
      </div>
      <div class="actions">
        <div class="button-group">
            <a class="minibutton tooltipped tooltipped-w"
               href="http://windows.github.com" aria-label="Open this file in GitHub for Windows">
                <span class="octicon octicon-device-desktop"></span> Open
            </a>
                <a class="minibutton js-update-url-with-hash"
                   href="/AugustHuang/MosaicHunter/edit/master/README.md"
                   data-method="post" rel="nofollow" data-hotkey="e">Edit</a>
          <a href="/AugustHuang/MosaicHunter/raw/master/README.md" class="minibutton " id="raw-url">Raw</a>
            <a href="/AugustHuang/MosaicHunter/blame/master/README.md" class="minibutton js-update-url-with-hash">Blame</a>
          <a href="/AugustHuang/MosaicHunter/commits/master/README.md" class="minibutton " rel="nofollow">History</a>
        </div><!-- /.button-group -->

            <a class="minibutton danger empty-icon"
               href="/AugustHuang/MosaicHunter/delete/master/README.md"
               data-method="post" data-test-id="delete-blob-file" rel="nofollow">

          Delete
        </a>
      </div><!-- /.actions -->
    </div>
      
  <div id="readme" class="blob instapaper_body">
    <article class="markdown-body entry-content" itemprop="mainContentOfPage"><h1>
<a name="user-content-mosaichunter" class="anchor" href="#mosaichunter" aria-hidden="true"><span class="octicon octicon-link"></span></a>MosaicHunter</h1>

<p>A script/tool for detecting postzygotic single-nucleotide mutations in human whole-genome sequencing data.</p>

<p>======Preparation</p>

<p>Make sure that you have installed the listed softwares, then added pre-installed softwares and the directory /your/MosaicHunter/directory/Tools in your PATH</p>

<p>In order to generate essential reference data and compile c, c++, java scripts, you should run this command once:</p>
<p>    seqpipe -m /your/MosaicHunter/directory/MosaicHunter.pipe preparation REFERENCE_DIR=/your/MosaicHunter/directory/Reference TOOLS_DIR=/your/MosaicHunter/directory/Tools</p>

<p>Pre-installed softwares required for the script:</p>
<p>    #bedtools: 2.15.0</p>
<p>    #samtools: 0.1.18</p>
<p>    #fastx_toolkit: 0.0.13</p>
<p>    #seqpipe: 0.4.12</p>
<p>    #blat</p>
<p>    #fastasplitn</p>

<p>Reference data for the script: (Please put them into /your/MosaicHunter/directory/Reference)</p>
<p>    #human_g1k_v37.fasta (available at <a href="http://soms.nibs.ac.cn:6237/human_g1k_v37.fasta">http://soms.nibs.ac.cn:6237/human_g1k_v37.fasta</a>)</p>
<p>    #human_g1k_v37.genome</p>
<p>    #human_hg19.fasta (available at <a href="http://soms.nibs.ac.cn:6239/human_hg19.fasta">http://soms.nibs.ac.cn:6239/human_hg19.fasta</a>)</p>
<p>    #all_repeats.b37.bed</p>
<p>    #PAR.b37.bed</p>
<p>    #dbsnp_137.b37.SNP_AF.tsv (available at <a href="http://soms.nibs.ac.cn:6235/dbsnp_137.b37.SNP_AF.tsv">http://soms.nibs.ac.cn:6235/dbsnp_137.b37.SNP_AF.tsv</a>)</p>
<p>    #observed_in_common.bed</p>

<p>Tools for the script: (Please put them into /your/MosaicHunter/directory/Tools)</p>
<p>    #generate_beta_log10_val_file.r</p>
<p>    #count_homopolymer.cpp</p>
<p>    #myjoin</p>
<p>    #my_join.pl</p>
<p>    #PileupFilter.java</p>
<p>    #genotyper.pipe</p>
<p>    #Yyx_genotype_log10lik_with_precalc_beta.c</p>
<p>    #Yyx_real_log10lik_from_baseQ.c</p>
<p>    #Yyx_individual_genotyper.c</p>
<p>    #LoFreq_call.c</p>
<p>    #sam2fa.pl</p>
<p>    #blat_best.pipe</p>
<p>    #highest-score.pl</p>
<p>    #calculate-score-coverage-identity.pl</p>
<p>    #intersect_bed12.pipe</p>
<p>    #my.grep</p>
<p>    #trimBamByBlock.pl</p>
<p>    #strand_bias.R</p>
<p>    #allele_pos_dist.R</p>
<p>    #splitSamByAllele.pl</p>

<p>======Run</p>

<p>To identify pSNM sites from the whole-genome sequencing data, you can run this command: </p>
<p>    seqpipe -m /your/MosaicHunter/directory/MosaicHunter.pipe MosaicHunter REFERENCE_DIR=/your/MosaicHunter/directory/Reference TOOLS_DIR=/your/MosaicHunter/directory/Tools TEMP_DIR=/your/temp/directory INPUT_BAM=example.bam INDEL_CNV_BED=example.bed PROJECT_NAME=example GENDER=M THREAD_NUM=5</p>
<p>        [INPUT_BAM]: the path of your input .bam file, the .bam file should be sorted and indexed</p>
<p>        [INDEL_CNV_BED]: the path of a .bed file containing all the CNV and indel-flanking(+-5bp) regions which will be masked in our pipeline</p>
<p>        [PROJECT_NAME]: a string used as the prefix and suffix of the output files's name</p>
<p>        [GENDER]:  the gender of the subject, F or M</p>
<p>        [THREAD_NUM]: the maximum number of threads for running the script</p>

<p>Recommended pre-processing of the .bam file:</p>
<p>    1) Removing the duplicated, improper-paired, and multi-hit reads</p>
<p>    2) Removing the reads with more than three mismatches</p>
<p>    3) Processing the reads by GATK's indel realignment and base quality score recalibration</p>


<p>To change the running order and the parameters of the Bayesian genotyepr and the error filters, you can edit the the scripts of MosaicHunter in /MosaicHunter/MosaicHunter.pipe, according to the user manual of seqpipe.</p>

<p>======Output</p>

<p>The final list of the pSNM candidates could be found at MosaicHunter_[PROJECT_NAME]/[PROJECT_NAME].mosaic.final.tsv</p>
<p>    The colunms in the final list represent: 1) chromosome</p>
<p>                                             2) position</p>
<p>                                             3) total depth</p>
<p>                                             4) reference nt</p>
<p>                                             5) alternative nt</p>
<p>                                             6) reference depth</p>
<p>                                             7) alternative depth</p>
<p>                                             8) -log10 of posterior probability of ref-hom genotype</p>
<p>                                             9) -log10 of posterior probability of het genotype</p>
<p>                                             10) -log10 of posterior probability of alt-hom genotype</p>
<p>                                             11) -log10 of posterior probability of mosaic genotype</p>
<p>                                             12) population allele fraction in dbSNP 137, -1 for annotated sites without information of allele fraction, -2 for unannotated sites</p>
<p>                                             13) sequence of +-500bp flanking regions</p></article>
  </div>

  </div>
</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" class="js-jump-to-line" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <form accept-charset="UTF-8" class="js-jump-to-line-form">
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" autofocus>
    <button type="submit" class="button">Go</button>
  </form>
</div>

        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div><!-- /.container -->
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer">
    <ul class="site-footer-links right">
      <li><a href="https://status.github.com/">Status</a></li>
      <li><a href="http://developer.github.com">API</a></li>
      <li><a href="http://training.github.com">Training</a></li>
      <li><a href="http://shop.github.com">Shop</a></li>
      <li><a href="/blog">Blog</a></li>
      <li><a href="/about">About</a></li>

    </ul>

    <a href="/">
      <span class="mega-octicon octicon-mark-github" title="GitHub"></span>
    </a>

    <ul class="site-footer-links">
      <li>&copy; 2014 <span title="0.05119s from github-fe130-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="/site/terms">Terms</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
        <li><a href="/contact">Contact</a></li>
    </ul>
  </div><!-- /.site-footer -->
</div><!-- /.container -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-fullscreen-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="fullscreen-contents js-fullscreen-contents" placeholder="" data-suggester="fullscreen_suggester"></textarea>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped tooltipped-w" aria-label="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped tooltipped-w"
      aria-label="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-x close js-ajax-error-dismiss" aria-label="Dismiss error"></a>
      Something went wrong with that request. Please try again.
    </div>


      <script crossorigin="anonymous" src="https://assets-cdn.github.com/assets/frameworks-df9e4beac80276ed3dfa56be0d97b536d0f5ee12.js" type="text/javascript"></script>
      <script async="async" crossorigin="anonymous" src="https://assets-cdn.github.com/assets/github-f56c7ffa6d9b791a74df3c50855fd5f54ed92203.js" type="text/javascript"></script>
      
      
        <script async src="https://www.google-analytics.com/analytics.js"></script>
  </body>
</html>

