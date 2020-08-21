// Get base URL
var loc = window.location;
var baseurl = loc.protocol + "//" + loc.host;

$(document).ready(function(){ 
    // Just to be sure
    $('body').css('display','none');
    // Fill in the included parts 
    includefiles();

    $('a[href^="#BASEURL#"]').attr('href',function(i,el){
        return relativetoabsolute(el);
    });

    // Make all img responsive
    $("img").each(function(index){$(this).addClass( "img-responsive") });
});

function includefiles(){
    var inLen = $('.include[file]').length;

    $('.include[file]').each(function (idx, el) {
        var included_file = $(this).attr('file');

        var wasinhead = $("head").html();
        $(this).load(included_file,function (response,status,xhr){
            if(status!=='success')
        {
            console.log('No sukcess');
        }
        if($(this).is("head"))
        {
            // console.log('This was head');
            $(this).html(function(i,origText){
                return  wasinhead + origText;
            });
            $(this).removeAttr('class');
        }

        $(this).find("img").attr('src',function(index,currVal){
            return relativetoabsolute(currVal);
        });

        $(this).find("a, link[href] script[href]").attr('href',function(index,currVal){
            return relativetoabsolute(currVal);
        });

        // Remove the include class and file attribute to avoid including it again
        $(this).removeAttr('file');
        $(this).removeClass('include');
        // This has to be inside load because it is called asynchronously
        // And we want to run it only after the last item was loaded
        if(idx==inLen-1)
        {  
            // Remove active item in menu
            $("nav .nav .active").removeClass("active");
            var chref = $(location).attr("href").substr(baseurl.length+1);

            var isTopLevelFile = chref.indexOf('.html') > -1 &&
                chref.indexOf('/') == -1 ;

            if(chref.length == 0){
                // Theat the home path differently
                $("nav .nav li").first().addClass("active"); 
            }
            else{
                if(isTopLevelFile){
                    $("nav .nav li a[href$='"+chref+"']")
                        .last().parent().addClass('active');
                }
                else
                {
                    // Get rid of the trailing filename
                    chref = chref.substr(chref,chref.lastIndexOf('/')).toLowerCase();
                    console.log("Current: " + chref);
                    
                    var array = ($("nav .nav li a:not([href$='.html'])").toArray());
                    var maxmatchel = array[0];
                    var maxmatch = -1;
                    var relpaths = $.map(array, function(val,ii){
                        var href = val.href.substr(baseurl.length+1).toLowerCase();
                        var tmpmatch = strcmp(href,chref);
                        if(tmpmatch>maxmatch){
                            maxmatch = tmpmatch;
                            maxmatchel = val;
                        }
                    });
                    $(maxmatchel).parent().addClass("active");

                   /* $("nav .nav li a:not([href$='.html'])").filter(function(){
                        var href = this.href.substr(baseurl.length+1);
                        if(href.length == 0){ return false;}
                        console.log("Menu:" + href);
                    return chref.indexOf(href,chref.length - href.length) !== -1;
                    }).last().parent().addClass("active");*/
                }
            }
            $('body').css({display:'block'});
        }
        });
    });

    if(inLen == 0)
    {
        // just make body visible and exit
        $('body').css({display:'block'});
    }
}

function relativetoabsolute(v){
    var n = v.lastIndexOf('#BASEURL#');
    if( n!= -1)
    {
        var valid = v.substring(n);
        var retval = valid.replace(/^#BASEURL#/,baseurl);

        //console.log('Old:' + v);
        //console.log('New:' + retval);
        return retval;
    }
}

function strcmp(a, b) {
    // retval is index of last matching character in strings
    // -1 means do not match at all
    a = a.toString(), b = b.toString();
    var retval = -1;
    for (;retval<Math.min(a.length, b.length);++retval){
        if(a.charAt(retval+1) !== b.charAt(retval+1)) return retval;
    }
    return retval;
}
