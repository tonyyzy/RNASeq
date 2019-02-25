
// SESSION DETAIL JAVASCRIPT
// $( document ).ready(function() {
//   var user_url = window.location.pathname
//   if (user_url.includes('session_detail')){ // executes if session_detail page loaded
//
// // var session_conditions = [];
// // for (var i = 0; i < 5; i++){
// //
// // }
// // $('#td_condition').text()
// });
// var conditions_max = [1,2,3,4,5,6,7,8,9,10];
// var conditions_count = [];
// for (var i = 0; i < conditions_max.length; i++){
//   if ($('#td_condition_[i]').text() != ''){ // checks if read_2 file upload is empty
//     console.log('td_condition_[i] passed');
//   }}
//
// else {
//   console.log('end')
// }
//
// $('#td_condition_[i]')

//
// $( ".td_condition" ).each(function( index ) {
//   return var myVar = $( this ).text();
// });
//
// $( ".td_replicate" ).each(function( index ) {
//   console.log( index + ": " + $( this ).text() );
// });


// SESSION CREATE JAVASCRIPT
$( document ).ready(function() {
  var user_url = window.location.pathname
  if (user_url.includes('session_create')){ // executes if session_create page loaded
    $('#div_id_organism').hide()
    $('#div_id_fasta_dna_file').hide()
    $('#div_id_fasta_cdna_file').hide()
    $('#div_id_gtf_file').hide()
    $('#div_id_accession').hide()
    $('#div_id_select_genome').hide()
    $('#div_id_genome').on('click', function(){
      // console.log('hi')
      if ($('#div_id_genome option:selected').val() == 'pre_index'){ // executes if pre_index selected from dropdown
        // console.log('great success with pre index')
        $('#div_id_organism').hide()
        $('#div_id_fasta_dna_file').hide()
        $('#div_id_fasta_cdna_file').hide()
        $('#div_id_gtf_file').hide()
        $('#div_id_accession').hide()
        $('#div_id_select_genome').show()
      }
      if ($('#div_id_genome option:selected').val() == 'user_provided'){ // executes if pre_index selected from dropdown
        // console.log('great success with pre index')
                $('#div_id_select_genome').hide()
        $('#div_id_organism').show()
        $('#div_id_fasta_dna_file').show()
        $('#div_id_fasta_cdna_file').show()
        $('#div_id_gtf_file').show()
        $('#div_id_accession').show()
      }
    })
  }
})




if ($('#div_id_libtype option:selected').val() == 'user_provided'){ // executes if user_provided selected from dropdown
  $('#div_id_read_1').show()
  $('#div_id_accession').show()
  $('#div_id_read_2').hide()
}







// SAMPLES CREATE JAVASCRIPT
$( document ).ready(function() {
  var user_url = window.location.pathname
  if (user_url.includes('samples_create')){ // executes if samples_create page loaded
    $('#div_id_read_1').hide()
    $('#div_id_read_2').hide()
    $('#div_id_accession').hide()
    $('#validate_library').hide()
    $('#div_id_libtype').on('click', function(){
      if ($('#div_id_libtype option:selected').val() == 'PE'){ // executes if PE selected from dropdown
        $('#div_id_read_1').show()
        $('#div_id_read_2').show()
        $('#div_id_accession').show()
      }
      if ($('#div_id_libtype option:selected').val() == 'SG'){ // executes if SG selected from dropdown
        $('#div_id_read_1').show()
        $('#div_id_accession').show()
        $('#div_id_read_2').hide()
      }
    })
  }
});

function validateLibrary(){
  if ($('#div_id_libtype option:selected').val() == 'PE'){ // executes if PE selected from dropdown
    if ($('#id_read_2').val() == ''){ // checks if read_2 file upload is empty
      $('#validate_library').show()
      return false;
    }
  }
}


// SAMPLES UPDATE JAVASCRIPT
$( document ).ready(function() {
  var user_url = window.location.pathname // stores current url as user_url
  if (user_url.includes('samples_update')){ // executes if samples_update page loaded
    $('#validate_library').hide()
    if ($('#div_id_libtype option:selected').val() == 'SG'){
    	$('#div_id_read_2').toggle()
    }
    $('#div_id_libtype').on('click', function(){
      if ($('#div_id_libtype option:selected').val() == 'PE'){ // executes if PE selected from dropdown
        $('#div_id_read_2').show()
        $('#read_2-clear_id').prop("checked", false);
      }
      if ($('#div_id_libtype option:selected').val() == 'SG'){ // executes if SG selected from dropdown
        $('#div_id_read_2').hide()
        $('#read_2-clear_id').prop("checked", true);
      }
    })
    }

});



// $('[name=options] option').filter(function() {
//     return ($(this).text() == 'Blue'); //To select Blue
// }).prop('selected', true);
//
//
// $('#div_id_libtype option:eq(2)')

// $('#div_id_libtype option:eq(2)').css("background-color", "yellow");
