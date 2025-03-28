<!DOCTYPE html>
<!--[if lt IE 7]><html class="no-js lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
<!--[if IE 7]><html class="no-js lt-ie9 lt-ie8"> <![endif]-->
<!--[if IE 8]><html class="no-js lt-ie9"> <![endif]-->
<!--[if gt IE 8]><!--><html lang="en" class="no-js"><!--<![endif]-->
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1">
    <title>Earthdata Login</title>
    <meta name="description" content="Earthdata Login">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">

    <!-- Google Tag Manager -->
    <script>(function(w,d,s,l,i){w[l]=w[l]||[];w[l].push(

      {'gtm.start': new Date().getTime(),event:'gtm.js'}

    );var f=d.getElementsByTagName(s)[0],
      j=d.createElement(s),dl=l!='dataLayer'?'&l='+l:'';j.async=true;j.src=
      'https://www.googletagmanager.com/gtm.js?id='+i+dl;f.parentNode.insertBefore(j,f);
    })(window,document,'script','dataLayer','GTM-WNP7MLF');</script>
    <!-- End Google Tag Manager -->

    <link href="https://cdn.earthdata.nasa.gov/eui/1.1.3/stylesheets/application.css" rel="stylesheet" />
    <link rel="stylesheet" href="/assets/application-0b92e5f8a20bb61e9b0121c043b3e00834a3fbb95cbe8a7bdabaf2a0dfcfc215.css" media="all" />
    <!--[if IE 7]>
      <link rel="stylesheet" href="/assets/font-awesome-ie7.min.css">
    <![endif]-->
    <link href="//netdna.bootstrapcdn.com/font-awesome/4.3.0/css/font-awesome.min.css" rel="stylesheet">
    <link href='https://fonts.googleapis.com/css?family=Source+Sans+Pro:300,700' rel='stylesheet' type='text/css'>
    <meta name="csrf-param" content="authenticity_token" />
<meta name="csrf-token" content="V6SA1LMdQa0q-eIFz851CqjlNhe1SBluk2mBrVx2sRepAKXX9iFdNlo5O74OIbaKss0h5d6QqFxXhR5B2_o2Jg" />
    

    <!-- Grid background: http://subtlepatterns.com/graphy/ -->
  </head>
  <body class="oauth authorize" data-turbolinks-eval=false>

    <!-- Google Tag Manager (noscript) -->
    <noscript>
      <iframe src="https://www.googletagmanager.com/ns.html?id=GTM-WNP7MLF"
                    height="0" width="0" style="display:none;visibility:hidden"></iframe>
    </noscript>
    <!-- End Google Tag Manager (noscript) -->

    <header id="earthdata-tophat2" style="height: 32px;"></header>
    <!--[if lt IE 7]>
      <p class="chromeframe">You are using an <strong>outdated</strong> browser. Please <a href="http://browsehappy.com/">upgrade your browser</a> or <a href="http://www.google.com/chromeframe/?redirect=true">activate Google Chrome Frame</a> to improve your experience.</p>
    <![endif]-->
    <header role="banner">
  <div class="container">
    <div id="masthead-logo">
      <h1><a class="ir" href="/">Earthdata Login</a></h1>
      <a href="/"><p class="masthead-text-logo">EARTHDATA <b>LOGIN</b></p></a>
    </div>
    <a id="hamburger" href="#"><img title="Mobile Menu" alt="Three horizontal lines stacked" src="/assets/hamburger-68c8505066427f3e3f6ee40b24cfd3c9f7c0fe93ee298b9046564637262115fa.png" /></a>
    <nav role="navigation" class="masthead">

      <div id="hide">
        <ul>
          <li><strong><a href="/documentation">Documentation</a></strong></li>
        </ul>
      </div>
    </nav>
  </div>
</header>

    <div class="container">
      







      <section id="callout-login">
  <div class="client-login">
    
    <br>
    <h3 class="client-description">
      Ocean Data Get File Authentication 
    </h3>

  </div>
  <form id="login" action="/login" accept-charset="UTF-8" method="post"><input name="utf8" type="hidden" value="&#x2713;" autocomplete="off" /><input type="hidden" name="authenticity_token" value="mftGDZpYSeAWi15vm-acFvASDdA7iATRmKRXDub2_uZnX2MO32RVe2ZLh9RaCV-W6joaIlBQteNcSMjiYXp51w" autocomplete="off" />
  <p>
    <label for="username">Username</label>
    <i class="fa fa-question-circle fa-question-circle--blue user-name" title="Login using either your Username or Email Address"></i>
    <input type="text" name="username" id="username" autofocus="autofocus" class="default" />
  </p>
  <p>
    <label for="password">Password</label><br />
    <input type="password" name="password" id="password" autocomplete="off" class="password" />
    <span class="password-showhide hide-password">
  </p>
  <p>
    <input type="hidden" name="client_id" id="client_id" value="pDPu0awH156XLrK6VV0Y0w" autocomplete="off" />
  </p>
  <p>
    <input type="hidden" name="redirect_uri" id="redirect_uri" value="https://oceandata.sci.gsfc.nasa.gov/getfile/urs/" autocomplete="off" />
  </p>
        <p><input type="hidden" name="response_type" id="response_type" value="code" autocomplete="off" /></p>
      <p><input type="hidden" name="state" id="state" autocomplete="off" /></p>
      <p><input type="checkbox" name="stay_in" id="stay_in" value="1" checked="checked" /> <label for="stay_in">Stay signed in (this is a private workstation)</label></p>

  <p class="button-with-notes">
    <input type="submit" name="commit" value="Log in" class="eui-btn--round eui-btn--green" data-disable-with="Log in" />
    <a class="eui-btn--round eui-btn--blue" id="register_new_account_button" href="/users/new?client_id=pDPu0awH156XLrK6VV0Y0w&amp;redirect_uri=https%3A%2F%2Foceandata.sci.gsfc.nasa.gov%2Fgetfile%2Furs%2F&amp;response_type=code">Register</a>
  </p>
  <p class="form-instructions">
    <em class="icon-question-sign"></em>
    <a class="" href="/retrieve_info">I don&rsquo;t remember my username</a>
    <br /><em class="icon-question-sign"></em>
    <a class="" href="/reset_passwords/new">I don&rsquo;t remember my password</a>
    <br />
    <em class="icon-question-sign"></em>
    <a href="javascript:feedback.showForm();" title = "Need Help? Click on the Feedback button to request help">Help</a>
  </p>

</form>
<aside class="govt-msg">
  <div class="nasa-logo"></div>
  <p>
    <strong>Why must I register?</strong>
  </p>
  <p>
    The Earthdata Login provides a single mechanism for user registration and profile management for all EOSDIS system components (DAACs, Tools, Services).
    Your Earthdata login also helps the EOSDIS program better understand the usage of EOSDIS services to improve user experience through customization of tools and improvement of services.
    EOSDIS data are openly available to all and free of charge except where governed by international agreements.
  </p>
</aside>

</section>
<section id="cta">
  <h3>Get single sign-on access to all your favorite EOSDIS sites</h3>
      <a class="eui-btn--round eui-btn--blue" href="/users/new?client_id=pDPu0awH156XLrK6VV0Y0w&amp;redirect_uri=https%3A%2F%2Foceandata.sci.gsfc.nasa.gov%2Fgetfile%2Furs%2F&amp;response_type=code">Register for a Profile</a>
</section>
<div class="govt-warning eui-info-box">
  <div class="warning-desktop">
    <p>
      <strong>
        Protection and maintenance of user profile information is described in
        <a href="https://www.nasa.gov/about/highlights/HP_Privacy.html">NASA's Web Privacy Policy.</a>
        </strong> 
    </p>
  </div>
  <div class="warning-mobile">
    <p>
      <strong>
        Protection and maintenance of user profile information is described in
            <a href="https://www.nasa.gov/about/highlights/HP_Privacy.html">NASA's Web Privacy Policy.</a>
      </strong> 
    </p>
  </div>
  <div class="warning-mobile-mini">
    <strong>
      US Govt Property. Unauthorized use subject to prosecution. Use subject to monitoring per
      <a href="https://nodis3.gsfc.nasa.gov/displayDir.cfm?t=NPD&c=2810&s=1E">NPD2810</a>.
    </strong>
  </div>
</div>

    </div>
    <footer role="contentinfo">
  <h3>For questions regarding the EOSDIS Earthdata Login, please contact <a href="javascript:feedback.showForm();" title="Earthdata Support form">Earthdata Support</a></h3>
  <ul>
    <li><b>V 4.218.0
</b></li>
    <li><a href="/">Home</a></li>
    <li><a href="/users/new">Register</a></li>
    <li><a title="NASA Home" href="http://www.nasa.gov">NASA</a></li>
    <li><a title="Accessability" href="https://www.nasa.gov/accessibility/">Accessibility</a></li>
  </ul>
  <p>NASA Official: Stephen Berrick</p>
</footer>

    <script src="/assets/application-d11b6707bbcb9bfd5966f378febb3c0d2ca25cae2f1a66b6f29b1b5178d9d045.js"></script>
    <script type="text/javascript">
    $(window).scroll(function(e){
      parallax();
    });
    function parallax(){
      var scrolled = $(window).scrollTop();
      $('#content').css('background-position', 'right ' + -(scrolled*0.25)+'px ');
    }
    </script>
    <script src="https://cdn.earthdata.nasa.gov/tophat2/tophat2.js" id="earthdata-tophat-script" data-show-fbm="true" data-show-status="true" data-status-api-url="https://status.earthdata.nasa.gov/api/v1/notifications" data-hide-daac-links="true"></script>
    <script type="text/javascript" src="https://fbm.earthdata.nasa.gov/for/URS4/feedback.js"></script>
    <script type="text/javascript">
      feedback.init();
    </script>
    <script type="text/javascript">
        setTimeout(function()
                {var a=document.createElement("script"); var b=document.getElementsByTagName("script")[0];
                    a.src=document.location.protocol+"//dnn506yrbagrg.cloudfront.net/pages/scripts/0013/2090.js?"+Math.floor(new Date().getTime()/3600000);
                    a.async=true;a.type="text/javascript";b.parentNode.insertBefore(a,b)}
                , 1);
    </script>

    <!-- BEGIN: DAP Google Analytics  -->
    <script language="javascript" id="_fed_an_ua_tag" src="https://dap.digitalgov.gov/Universal-Federated-Analytics-Min.js?agency=NASA&subagency=GSFC&dclink=true"></script>
    <!-- END: DAP Google Analytics  -->

    
  </body>
</html>
