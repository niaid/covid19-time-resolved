
// Uses Github API to retrieve a list of repositories by the lab

var githubAPI = "https://api.github.com/search/repositories?q=";
var repositories = [];
var outhtml = "";
$(function() {
	$.when(repoSearch("org:niaid topic:niaid-tsang-lab"),repoSearch("org:niaid protocol"),repoSearch("baseline user:kotliary"),
	repoSearch("topic:niaid-tsang-lab user:MattPM")).done(function(a1,a2,a3) {
		outputPageContent();
	});
});

function repoSearch(query) {
    return $.ajax({
    	dataType: "json",
    	url: githubAPI + query,
	success: function(data) {
	    if (data.total_count > 0)
	    	repositories = repositories.concat(data.items)
	},
   	headers: {'Accept':'application/vnd.github.mercy-preview+json'}
    })
} // end repoSearch()

function outputPageContent() {
     if(repositories.length == 0) { outhtml = outhtml + '<p>No repositories found!</p>'; }
     else {
         $.each(repositories.sort(compare), function(index,repo) {
	     var name = repo.name.replace(/_/g," ");
	     var desc = repo.description;
 	     var url = repo.html_url;
	     if (desc == null)
		desc = "Description Not Available";
             outhtml = outhtml + '<div class="col-xs-12 col-sm-6 col-md-3"><div class="thumbnail block-codes ">' + 
             	'<h4 class="codes-title"><a class="group-link" href="' + url + '" target="blank">' + 
		name + '</a></h4><hr><p class="codes-desc">' + desc + '</p>' +
            	'<div class="codes-bottom"><div class="codes-bottom-row"><div class="col-xs-12">' +
                '<p><a href="' + url + '" class="btn btn-primary codes-download btn-block " target="blank" role="button">Go to repository</a></p>' +
                '</div> <!-- /col-xs-12 --></div><!-- /codes-bottom-row--></div> <!-- /codes-bottom --></div></div><!-- /col-xs-12 col-sm-6 col-md-3 -->'   
         });
     }
     $('#dynamic-list').html(outhtml);
} // end outputPageContent()

function compare(a, b) {
  // Use toUpperCase() to ignore character casing
  const bandA = a.name.toUpperCase();
  const bandB = b.name.toUpperCase();

  let comparison = 0;
  if (bandA > bandB) {
    comparison = 1;
  } else if (bandA < bandB) {
    comparison = -1;
  }
  return comparison;
} // end compare
    