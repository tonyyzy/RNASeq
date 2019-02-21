
// $('[name=options] option').filter(function() {
//     return ($(this).text() == 'Blue'); //To select Blue
// }).prop('selected', true);
//
//
// $('#div_id_libtype option:eq(2)')

// $('#div_id_libtype option:eq(2)').css("background-color", "yellow");



// SAMPLES CREATE JAVASCRIPT
$( document ).ready(function() {
  var user_url = window.location.pathname
  if (user_url.includes('samples_create')){ // executes if samples_create page loaded
    $('#div_id_read_1').hide()
    $('#div_id_read_2').hide()
    $('#div_id_accession').hide()
    $('#div_id_libtype').on('click', function(){
      if ($('#div_id_libtype option:selected').val() == 'PE'){ // executes if PE selected from dropdown
        $('#div_id_read_1').show()
        $('#div_id_read_2').show()
        $('#div_id_accession').show()
      }
      if ($('#div_id_libtype option:selected').val() == 'SG'){ // executes if SG selected from dropdown
        $('#div_id_read_2').hide()
      }
    })
  }
  else {
    console.log('no blinding required')
  }
});

// SAMPLES UPDATE JAVASCRIPT
$( document ).ready(function() {
  var user_url = window.location.pathname // stores current url as user_url
  if (user_url.includes('samples_update')){ // executes if samples_update page loaded
    if ($('#div_id_libtype option:selected').val() == 'SG'){
    	$('#div_id_read_2').hide()
    }
    $('#div_id_libtype').on('click', function(){
      if ($('#div_id_libtype option:selected').val() == 'PE'){ // executes if PE selected from dropdown
        $('#div_id_read_2').show()
        $('#read_2-clear_id').prop("checked", false);
      }
      if ($('#div_id_libtype option:selected').val() == 'SG'){ // executes if PE selected from dropdown
        $('#div_id_read_2').hide()
        $('#read_2-clear_id').prop("checked", true);
      }
    })
    }
  else {
    console.log('no hiding required')
  }
});


function validateLibrary(){
  if ($('#div_id_libtype option:selected').val() == 'PE'){ // executes if PE selected from dropdown
    if ($('#id_read_2').val() == ''){
      alert('please select second file for paired end read submission');
      return false;
    }
  }
}
