#/*@@
#  @file      CreateImplementationBindings.pl
#  @date      Sun Jul  4 17:09:54 1999
#  @author    Tom Goodale
#  @desc
#
#  @enddesc
#@@*/

sub CreateImplementationBindings
{
  my($bindings_dir, $rhparameter_db, $rhinterface_db, $configuration_db) = @_;
  my($i, $start_dir, $thorn);
  my(@data, @thorns, @ancestors, @friends, @requires_thorns, @activates_thorns);

  if(! $build_dir)
  {
    $build_dir = "$bindings_dir/build";
  }

  if(! -d $bindings_dir)
  {
    mkdir("$bindings_dir", 0755) || die "Unable to create $bindings_dir";
  }
  $start_dir = `pwd`;

  chdir $bindings_dir;

  if(! -d 'Implementations')
  {
    mkdir('Implementations', 0755) || die 'Unable to create Implementations directory';
  }

  if(! -d 'include')
  {
    mkdir('include', 0755) || die 'Unable to create include directory';
  }

  @thorns = sort split(' ', $rhinterface_db->{'THORNS'});
  if(! -d 'include')
  {
    mkdir('include', 0755) || die 'Unable to create include directory';
  }
  if(! -d 'include/CactusBindings')
  {
    mkdir('include/CactusBindings', 0755) || die 'Unable to create include/CactusBindings directory';
  }
  my $thorn;
  foreach $thorn (@thorns) {
    if(! -d "include/$thorns") {
      mkdir("include/$thorn",0755) || die "Unable to create include/$thorn directory";
    }
  }

  @data = map { "void CCTKi_BindingsThorn_$_(void);" } @thorns;
  push(@data, '');

  push(@data, 'int CCTKi_BindingsImplementationsInitialise(void);');
  push(@data, 'int CCTKi_BindingsImplementationsInitialise(void)');
  push(@data, '{');
  push(@data, map { "  CCTKi_BindingsThorn_$_();" } @thorns);
  push(@data, '  return 0;');
  push(@data, '}');
  push(@data, "\n");  # workaround for perl 5.004_04 to add a trailing newline

  $dataout = join ("\n", @data);
  &WriteFile('Implementations/ImplementationBindings.c',\$dataout);

  $dataout = 'SRCS = ImplementationBindings.c';
  &WriteFile('Implementations/make.code.defn',\$dataout);

  if(! -d "$build_dir")
  {
    mkdir("$build_dir", 0755) || die "Unable to create $build_dir";
  }

  chdir "$build_dir";

  foreach $thorn (@thorns)
  {
    if(! -d "$thorn")
    {
      mkdir("$thorn", 0755) || die "Unable to create $build_dir/$thorn";
    }

    $myimp = $rhinterface_db->{"\U$thorn\E IMPLEMENTS"};
    @ancestors = map { "    \"$_\"," } split (' ', $rhinterface_db->{"IMPLEMENTATION \U$myimp\E ANCESTORS"});
    @friends = map { "    \"$_\"," } split (' ', $rhinterface_db->{"\U$thorn\E FRIEND"});
    @requires_thorns = map { "    \"$_\"," } split (' ', $configuration_db->{"\U$thorn\E REQUIRES THORNS"});
    @activates_thorns = map { "    \"$_\"," } split (' ', $configuration_db->{"\U$thorn\E ACTIVATES THORNS"});

    @data = ();
    push(@data, '#include <stdio.h>');
    push(@data, '');
    push(@data, '#include "cctki_ActiveThorns.h"');
    push(@data, '');

    push(@data, "void CCTKi_BindingsThorn_${thorn}(void);");
    push(@data, "void CCTKi_BindingsThorn_${thorn}(void)");
    push(@data, '{');
    push(@data, "  const char *name[] = {\"$thorn\", 0};");
    push(@data, "  const char *implementation[] = {\"$myimp\", 0};");
    $i = 3;
    if (@ancestors)
    {
      push(@data, '  const char *ancestors[] =');
      push(@data, '  {');
      push(@data, @ancestors);
      push(@data, '    0,');
      push(@data, '  };');
      push(@data, '');
      $i++;
    }

    if (@friends)
    {
      # Just pass the ones this thorn has declared itself to be friends with.
      push(@data, '  const char *friends[] =');
      push(@data, '  {');
      push(@data, @friends);
      push(@data, '    0,');
      push(@data, '  };');
      push(@data, '');
      $i++;
    }

    if (@requires_thorns)
    {
      push(@data, '  const char *requires_thorns[] =');
      push(@data, '  {');
      push(@data, @requires_thorns);
      push(@data, '    0,');
      push(@data, '  };');
      push(@data, '');
      $i++;
    }

    if (@activates_thorns)
    {
      push(@data, '  const char *activates_thorns[] =');
      push(@data, '  {');
      push(@data, @activates_thorns);
      push(@data, '    0,');
      push(@data, '  };');
      push(@data, '');
      $i++;
    }

    push(@data, '  /*');
    push(@data, '   * Should be able to do below with a constant initializer');
    push(@data, '   * but sr8000 compiler doesn\'t like it.');
    push(@data, '   * So have to laboriously assign values to each member of array.');
    push(@data, '   */');
    push(@data, "  struct iAttributeList attributes[$i];");
    push(@data, '');
    push(@data, '  attributes[0].attribute                = "name";');
    push(@data, '  attributes[0].AttributeData.StringList = name;');
    push(@data, '  attributes[1].attribute                = "implementation";');
    push(@data, '  attributes[1].AttributeData.StringList = implementation;');
    $i = 2;
    if (@ancestors)
    {
      push(@data, "  attributes[$i].attribute                = \"ancestors\";");
      push(@data, "  attributes[$i].AttributeData.StringList = ancestors;");
      $i++;
    }
    if (@friends)
    {
      push(@data, "  attributes[$i].attribute                = \"friends\";");
      push(@data, "  attributes[$i].AttributeData.StringList = friends;");
      $i++;
    }
    if (@requires_thorns)
    {
      push(@data, "  attributes[$i].attribute                = \"requires thorns\";");
      push(@data, "  attributes[$i].AttributeData.StringList = requires_thorns;");
      $i++;
    }
    if (@activates_thorns)
    {
      push(@data, "  attributes[$i].attribute                = \"activates thorns\";");
      push(@data, "  attributes[$i].AttributeData.StringList = activates_thorns;");
      $i++;
    }
    push(@data, "  attributes[$i].attribute                = 0;");
    push(@data, "  attributes[$i].AttributeData.StringList = 0;");
    push(@data, '');
    push(@data, '  CCTKi_RegisterThorn(attributes);');
    push(@data, '}');
    push(@data, "\n");  # workaround for perl 5.004_04 to add a trailing newline

    $dataout = join ("\n", @data);
    &WriteFile("$thorn/cctk_ThornBindings.c",\$dataout);

    $dataout = 'SRCS = cctk_ThornBindings.c';
    &WriteFile("$thorn/make.code.defn",\$dataout);
  }

  chdir($start_dir);
}

1;
